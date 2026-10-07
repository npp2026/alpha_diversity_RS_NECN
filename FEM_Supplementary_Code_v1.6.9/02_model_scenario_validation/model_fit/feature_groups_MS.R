# Grouping inherited from v1.5.2; selection itself is inside training folds.
define_feature_groups <- function(data) {
  
  cols <- names(data)
  
  # === Define ID/metadata and target columns (to exclude) ===
  id_cols <- c("plot_id", "plot_year", "aug_offset", "image_year", "year", 
               "Lat_Export", "Lon_Export", "site", "system:index", ".geo")
  target_cols <- c("Rich_tree", "Shannon_wiener", "Simpson",
                   "Herb_rich", "Herb_shannon", "Herb_simpson", "Herb_pielou",
                   "Shrub_rich", "Shrub_shannon", "Shrub_simpson", "Shrub_pielou",
                   "Shrub_density", "Treelet_density")
  
  # === Plot/Site-level attributes (PL_site) ===
  # These are field-measured or pre-defined site characteristics
  pl_site_candidates <- c("Forest_age", "Forest_type", "PFT", "AGB")
  group_PL_site <- intersect(pl_site_candidates, cols)
  
  # === Human Footprint (Human_footprint) ===
  # Human influence/disturbance indices
  hf_candidates <- c("HPI")
  group_Human_footprint <- intersect(hf_candidates, cols)
  
  # =========================================================================
  # === Productivity (Prod) ===
  # 生产力假说: 能量/资源可用性 → 物种数
  # 包含: 光谱值、植被指数、MODIS生产力指标 (均值/最大值/总和)
  # =========================================================================
  
  # Landsat bands & VIs (instantaneous seasonal values)
  raw_patterns <- c("^B[1-3]_", "^NIR_", "^SWIR[12]_",
                    "^NDVI_", "^EVI_", "^SAVI_", "^NBR_", "^NDMI_", "^GRVI_", "^VARI_", "^NDSI_")
  vi_instant <- grep(paste(raw_patterns, collapse = "|"), cols, value = TRUE)
  
  exclude_patterns <- "_(stdDev|CV|contrast|ent|asm|corr|var|diss|SDil|SEro|SGrad|TSD|GS|Diff|Amp|Ent|Range)"
  prod_instant <- vi_instant[!grepl(exclude_patterns, vi_instant)]
  
  # Landsat STMs magnitude (时间序列的均值/中位数/极值 - 代表典型水平)
  stms_magnitude <- grep("_(Mean|Med|P10|P90|Max|Min)_GS$", cols, value = TRUE)
  
  # MODIS Productivity magnitude (GPP, LAI, fPAR mean/max/sum)
  modis_prod_magnitude <- grep("^Prod_(GPP|LAI|fPAR)_(Mean|Max|Sum)$", cols, value = TRUE)
  
  group_Prod <- unique(c(prod_instant, stms_magnitude, modis_prod_magnitude))
  
  # =========================================================================
  # === Heterogeneity (Het) - 纯空间维度 ===
  # 空间异质性假说 (SVH): 空间复杂度 → 微生境分化 → 生态位分化
  # GEE计算: reduceNeighborhood() - 空间邻域统计
  # 包含: 空间CV/SD、纹理、形态学、MODIS空间SD
  # =========================================================================
  
  # Landsat 空间统计 (邻域窗口CV/SD)
  # GEE: image.reduceNeighborhood({kernel: ee.Kernel.square(radius)})
  het_spatial_stats <- grep("_(stdDev|CV)_w[0-9]+|_(stdDev|CV)_[0-9]+m|^SG_|^SD_|^aSG_", cols, value = TRUE)
  
  # Landsat 纹理特征 (GLCM)
  # GEE: image.glcmTexture() → 空间纹理
  het_texture <- grep("_(contrast|ent|asm|corr|var|diss)_|_Ent_", cols, value = TRUE)
  
  # Landsat 形态学特征 (膨胀/腐蚀/梯度)
  het_morph <- grep("_(SDil|SEro|SGrad)_", cols, value = TRUE)
  
  # MODIS 物候空间异质性 (邻域SD/Range/Diversity)
  # GEE: phenoStack.reduceNeighborhood({kernel: kernel_1km/2km})
  # 保持宽松pattern以兼容未来可能的新尺度
  pheno_spatial_het <- grep("^Pheno_.*(SD_|Range_|Diversity)", cols, value = TRUE)
  
  # MODIS 生产力空间异质性 (邻域SD/Range) [v3.8修正: 移除CV]
  # GEE: prodStack.reduceNeighborhood({kernel: kernel_1km/2km})
  modis_prod_spatial <- grep("^Prod_(GPP|LAI|fPAR)_(SD_[12]km|Range_[12]km)$", cols, value = TRUE)
  
  group_Het <- unique(c(het_spatial_stats, het_texture, het_morph, 
                        pheno_spatial_het, modis_prod_spatial))
  
  # =========================================================================
  # === Temporal (Temp) - 时间维度 ===
  # 时间动态假说: 时间变异 → 竞争重置/资源分化 → 物种共存
  # GEE计算: collection.reduce(stdDev/CV) - 时间序列统计
  # 包含: 时间CV/SD、季节差异、物候参数、MODIS时间CV
  # =========================================================================
  
  # Landsat 时间序列统计 (目标年 ±1 年影像的时间stdDev)
  # GEE: collection.reduce(ee.Reducer.stdDev()) → _Std_GS, _TSD_
  time_ts_stats <- grep("_(Std|Range)_GS$|_TSD_", cols, value = TRUE)
  
  # Landsat 季节差异 (物候转换强度)
  # GEE: summerPeak.subtract(springPeak) → SumSprDiff, SeasonAmp
  time_seasonal_diff <- grep("(SumSprDiff|AutSumDiff|SeasonAmp)$", cols, value = TRUE)
  
  # MODIS 物候参数 (时间模式/季节特征)
  # GEE: phenology timing, duration, amplitude
  pheno_timing <- grep("^Pheno_(Greenup|Peak|Senescence|GSL|GreenupDur|SenescenceDur|Asymmetry|Amp)$", cols, value = TRUE)
  
  # MODIS 生产力时间变异 (3年时间序列CV) [v3.8修正: 从Het移入]
  # GEE: fpar_lai_col.reduce(ee.Reducer.stdDev()) / mean → CV
  modis_prod_temporal <- grep("^Prod_(GPP|LAI|fPAR)_CV$", cols, value = TRUE)
  
  group_Temp <- unique(c(time_ts_stats, time_seasonal_diff, pheno_timing, modis_prod_temporal))
  
  # =========================================================================
  # === Environment (Env) ===
  # 环境决定论: 非生物环境 → 物种分布
  # 包含: 地形、气候、土壤
  # =========================================================================
  
  # Topography core
  topo_core <- intersect(c("ELEV", "Elev", "Elevation", "Slope", "Northness", "Eastness", 
                           "TPI_1km", "TRI", "Aspect", "Topographic_Wetness_Index"), cols)
  # Topography heterogeneity
  topo_het <- grep("^(Elev|Slope|TPI_1km)_CV_", cols, value = TRUE)
  group_Topo <- unique(c(topo_core, topo_het))
  
  # Climate - Updated to include Seasonality_ patterns
  clim_patterns <- "^BIO[0-9]+(_|$)|^SPEI_|^CL_|^Clim_|^Base_|_T$|_P$|^Annual_|^Coldest_|^Warmest_|^Driest_|^Wettest_|^Seasonality_|^Isothermality_|^Diurnal_"
  group_Clim <- grep(clim_patterns, cols, value = TRUE)
  
  # Soil
  soil_patterns <- c("^BD$", "^SOC$", "^STN$", "^STP$", "^STK$", "^pH$",
                     "^CEC$", "^CF$", "^Btcly$", "^Btslt$", "^Btsnd$", "^Texcls$", "^Thickness$")
  group_Soil <- grep(paste(soil_patterns, collapse = "|"), cols, value = TRUE)
  
  group_Env <- unique(c(group_Topo, group_Clim, group_Soil))
  
  # v1.6.1: source-removal requires disjoint sources. Seasonal contrasts
  # were also matched by the broad productivity prefix; terrain CV matched Het.
  group_Het <- setdiff(group_Het, group_Env)
  group_Prod <- setdiff(group_Prod, c(group_Env, group_Het, group_Temp))

  # =========================================================================
  # === Assemble return object ===
  # =========================================================================
  groups <- list(
    # ===== Main groups (6 groups) =====
    Env = group_Env,
    Prod = group_Prod,
    Het = group_Het,      # v3.8: 纯空间维度
    Temp = group_Temp,    # v3.8: 时间维度 + 物候
    PL_site = group_PL_site,
    Human_footprint = group_Human_footprint,
    
    # ===== Environment subgroups =====
    Env_Topo = group_Topo,
    Env_Clim = group_Clim,
    Env_Soil = group_Soil,
    
    # ===== Productivity subgroups =====
    Prod_Landsat = unique(c(prod_instant, stms_magnitude)),
    Prod_MODIS = modis_prod_magnitude,
    
    # ===== Heterogeneity subgroups (空间维度) =====
    Het_Spatial_Stats = het_spatial_stats,   # 空间CV/SD
    Het_Texture = het_texture,               # 纹理
    Het_Morph = het_morph,                   # 形态学
    Het_Pheno_Spatial = pheno_spatial_het,   # 物候空间SD
    Het_MODIS_Spatial = modis_prod_spatial,  # 生产力空间SD [v3.8重命名]
    
    # ===== Temporal subgroups (时间维度) =====
    Temp_TS_Stats = time_ts_stats,           # 时间序列统计
    Temp_Seasonal = time_seasonal_diff,      # 季节差异
    Temp_Pheno = pheno_timing,               # 物候参数
    Temp_MODIS_CV = modis_prod_temporal,     # 生产力时间CV [v3.8新增]
    
    # ===== Legacy aliases (兼容旧代码) =====
    Het_Stats = het_spatial_stats,           # 别名
    Het_Pheno = pheno_spatial_het,           # 别名
    Het_MODIS_Prod = modis_prod_spatial,     # 别名 (注意: 不再包含CV)
    Temp_Landsat = unique(c(time_ts_stats, time_seasonal_diff)),
    Temp_MODIS = unique(c(pheno_timing, modis_prod_temporal))  # v3.8: 包含时间CV
  )
  
  # === Check for ungrouped ===
  all_grouped <- unique(unlist(groups[c("Env", "Prod", "Het", "Temp", "PL_site", "Human_footprint")]))
  feature_cols <- setdiff(cols, c(id_cols, target_cols))
  ungrouped <- setdiff(feature_cols, all_grouped)
  
  # Print statistics
  cat("\n", strrep("=", 70), "\n")
  cat(" 特征分组统计 (v3.8 - 时空维度精确分类)\n")
  cat(strrep("=", 70), "\n")
  for (g in c("Env", "Prod", "Het", "Temp", "PL_site", "Human_footprint")) {
    cat(sprintf("  %-16s: %4d 特征\n", g, length(groups[[g]])))
  }
  cat(strrep("-", 70), "\n")
  cat(" 子组详情:\n")
  cat(sprintf("  Prod_Landsat:      %4d | Prod_MODIS:        %4d\n", 
              length(groups$Prod_Landsat), length(groups$Prod_MODIS)))
  cat(strrep("-", 70), "\n")
  cat(" Het子组 (空间维度 - reduceNeighborhood):\n")
  cat(sprintf("  Het_Spatial_Stats: %4d | Het_Texture:       %4d\n", 
              length(groups$Het_Spatial_Stats), length(groups$Het_Texture)))
  cat(sprintf("  Het_Morph:         %4d | Het_Pheno_Spatial: %4d | Het_MODIS_Spatial: %4d\n", 
              length(groups$Het_Morph), length(groups$Het_Pheno_Spatial), length(groups$Het_MODIS_Spatial)))
  cat(strrep("-", 70), "\n")
  cat(" Temp子组 (时间维度 - collection.reduce + 物候):\n")
  cat(sprintf("  Temp_TS_Stats:     %4d | Temp_Seasonal:     %4d\n", 
              length(groups$Temp_TS_Stats), length(groups$Temp_Seasonal)))
  cat(sprintf("  Temp_Pheno:        %4d | Temp_MODIS_CV:     %4d\n", 
              length(groups$Temp_Pheno), length(groups$Temp_MODIS_CV)))
  cat(strrep("-", 70), "\n")
  cat(" Env子组:\n")
  cat(sprintf("  Env_Topo:          %4d | Env_Clim:          %4d | Env_Soil: %4d\n", 
              length(groups$Env_Topo), length(groups$Env_Clim), length(groups$Env_Soil)))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  %-16s: %4d 特征 (六大组去重后)\n", "总计", length(all_grouped)))
  
  # Report ungrouped
  if (length(ungrouped) > 0) {
    cat(strrep("=", 70), "\n")
    cat(sprintf(" 警告: %d 个特征未分组:\n", length(ungrouped)))
    cat("  ", paste(ungrouped, collapse = ", "), "\n")
  } else {
    cat(sprintf(" ✓ 所有 %d 个特征已分组\n", length(feature_cols)))
  }
  cat(strrep("=", 70), "\n")
  
  # v3.8: 添加分组维度说明
  cat("\n [v3.8分组说明]\n")
  cat("  Het = 纯空间维度 (reduceNeighborhood): 空间CV/SD, 纹理, MODIS空间SD\n")
  cat("  Temp = 时间维度 (collection.reduce): 时间CV/SD, 季节差异, 物候, MODIS时间CV\n")
  cat("  注意: Prod_*_CV 已从Het移入Temp (时间序列CV, 非空间CV)\n\n")
  
  return(groups)
}
