/**
 * Code S2. MODIS productivity and phenology extraction.
 *
 * Purpose
 *   Extract MODIS fPAR, LAI, GPP, phenology, and MODIS-derived spatial
 *   heterogeneity metrics for model training.
 *
 * Method alignment
 *   - Target-year +/- 1 three-year temporal window.
 *   - Productivity: MOD15A2H fPAR/LAI and MOD17A2H GPP, May-September.
 *   - Phenology: MCD12Q2 transition timing and EVI amplitude.
 *   - Product-specific QA masks, scale factors, and valid-range filters.
 *   - Phenological transition dates converted from epoch days to day of year.
 *   - MODIS values are sampled onto the 30 m analysis grid with default
 *     nearest-neighbour assignment, preserving the source-product support.
 */

// ===================================================================================
// 1. Configuration
// ===================================================================================
var CONFIG = {
  assetId: '',
  
  augmentation: {
    enabled: false,
    offsets: [0]
  },
  
  // Three-year window: target year +/- 1 year.
  yearBufferForComposite: 1,
  
  modis: {
    scales: {
      small: 1000,  // 1km
      large: 2000   // 2km
    }
  },
  
  productivity: {
    enabled: true,
    startMonth: 5,   // May
    endMonth: 9,     // September
    fparMax: 1.0,    // Valid fPAR range: 0-1
    laiMax: 10.0,    // Valid LAI range: 0-10 m²/m²
    gppMax: 0.1      // Valid GPP range: 0-0.1 kg*C/m² per 8-day
  },
  
  export: { 
    folder: 'plot_predictor_exports', 
    scale: 30 
  },
  
  tileScale: 4
};

var TEST_MODE = false;
var ALL_YEARS = TEST_MODE ? [2010] : [2008, 2009, 2010, 2011, 2012, 2014, 2015, 2016, 2017];

// ===================================================================================
// 2. Setup
// ===================================================================================
var albers_wkt = 'PROJCS["Asia_North_Albers_Equal_Area_Conic", ' +
  'GEOGCS["GCS_WGS_1984", DATUM["D_WGS_1984", SPHEROID["WGS_1984",6378137,298.257223563]], ' +
  'PRIMEM["Greenwich",0], UNIT["Degree",0.0174532925199433]], ' +
  'PROJECTION["Albers"], PARAMETER["False_Easting",0], PARAMETER["False_Northing",0], ' +
  'PARAMETER["Central_Meridian",105], PARAMETER["Standard_Parallel_1",25], ' +
  'PARAMETER["Standard_Parallel_2",47], PARAMETER["Latitude_Of_Origin",30], ' +
  'UNIT["Meter",1]]';

var plots = ee.FeatureCollection(CONFIG.assetId).map(function(f) {
  var lat = ee.Number(f.get('lat'));
  var lon = ee.Number(f.get('long'));
  return f.setGeometry(ee.Geometry.Point([lon, lat]))
          .set('Lat_Export', lat)
          .set('Lon_Export', lon);
});

// ===================================================================================
// 3. MODIS QA/QC Masking Helpers
// ===================================================================================
function maskMCD12Q2(image) {
  // QA_Overall_1: keep best/good main-cycle retrievals, then drop common fill dates.
  var qaOverall = image.select('QA_Overall_1');
  var qaMask = qaOverall.lte(1);

  var dateBands = image.select(['Greenup_1', 'Peak_1', 'Senescence_1', 'Dormancy_1']);
  var validDates = dateBands.reduce(ee.Reducer.min()).gt(0)
    .and(dateBands.reduce(ee.Reducer.max()).lt(32766));

  return image.updateMask(qaMask.and(validDates));
}

function maskMOD15A2H(image) {
  // FparLai_QC bit 0 is MODLAND_QC; 0 indicates good quality.
  var qc = image.select('FparLai_QC');
  var good = qc.bitwiseAnd(1).eq(0);
  return image.updateMask(good);
}

function maskMOD17A2H(image) {
  // Psn_QC bit 0 is MODLAND_QC; 0 indicates good quality.
  var qc = image.select('Psn_QC');
  var good = qc.bitwiseAnd(1).eq(0);
  return image.updateMask(good);
}

// ===================================================================================
// 4. Helper Function: Convert Epoch Days to DOY
// ===================================================================================
/**
 * Convert MODIS phenology dates (days since 1970-01-01) to Day of Year (DOY)
 * 
 * MODIS MCD12Q2 stores transition dates as days since Unix epoch (1970-01-01).
 * This function converts to DOY for the target year.
 * 
 * @param {ee.Image} epochDaysImage - Image with values as days since 1970-01-01
 * @param {ee.Number} targetYear - The target year for DOY calculation
 * @returns {ee.Image} - Image with DOY values (1-365 or 1-366 for leap years)
 * 
 * Formula: DOY = epochDays - epochDaysOfJan1 + 1
 * Example: If epochDays = 18628 and Jan 1, 2021 = 18628, then DOY = 1
 */
function epochDaysToDOY(epochDaysImage, targetYear) {
  // Days from 1970-01-01 to Jan 1 of target year
  var jan1 = ee.Date.fromYMD(targetYear, 1, 1);
  var epochOffset = jan1.difference(ee.Date('1970-01-01'), 'day');
  
  // DOY = epochDays - epochOffset + 1.
  return epochDaysImage.subtract(epochOffset).add(1);
}

function prepareAnnualPhenology(img) {
var ownYear = ee.Date(img.get('system:time_start')).get('year');
var g = epochDaysToDOY(img.select('Greenup_1'), ownYear);
var p = epochDaysToDOY(img.select('Peak_1'), ownYear);
var s = epochDaysToDOY(img.select('Senescence_1'), ownYear);
var d = epochDaysToDOY(img.select('Dormancy_1'), ownYear);
var amplitudeRaw = img.select('EVI_Amplitude_1');
var amplitude = amplitudeRaw.multiply(0.0001)
  .updateMask(amplitudeRaw.gte(0).and(amplitudeRaw.lte(10000)));
return ee.Image.cat([g, p, s, d.subtract(g), p.subtract(g),
  d.subtract(s), p.subtract(g).subtract(d.subtract(s)), amplitude])
  .rename(['greenup', 'peak', 'senescence', 'gsl', 'greenupDuration',
    'senescenceDuration', 'asymmetry', 'amp']);
}

// ===================================================================================
// 5. MODIS Phenology Function (17 features)
// ===================================================================================
function getModisPhenology(year) {
  var outputBandNames = [
    'Pheno_Greenup', 'Pheno_Peak', 'Pheno_Senescence', 'Pheno_Amp',
    'Pheno_GSL', 'Pheno_GreenupDur', 'Pheno_SenescenceDur', 'Pheno_Asymmetry',
    'Pheno_Greenup_SD_1km', 'Pheno_GSL_SD_1km', 'Pheno_Peak_SD_1km', 
    'Pheno_Amp_SD_1km', 'Pheno_Greenup_Range_1km',
    'Pheno_Greenup_SD_2km', 'Pheno_GSL_SD_2km',
    'Pheno_Diversity_1km', 'Pheno_Diversity_2km'
  ];
  
  // Three-year window: target year +/- 1 year.
  var yearNum = ee.Number(year);
  var startYear = yearNum.subtract(CONFIG.yearBufferForComposite);
  var endYear = yearNum.add(CONFIG.yearBufferForComposite);
  
  var start = ee.Date.parse('YYYY-MM-dd', startYear.format('%d').cat('-01-01'));
  var end = ee.Date.fromYMD(endYear.add(1), 1, 1); // exclusive end
  
  var modisCol = ee.ImageCollection('MODIS/061/MCD12Q2')
    .filterDate(start, end)
    .map(maskMCD12Q2);
  
  var colSize = modisCol.size();
    
  return ee.Algorithms.If(
    colSize.gt(0),
    (function() {
      // Convert every annual image relative to its OWN product year before
      // averaging. Averaging epoch days first is biased when a year is masked.
      var annual = modisCol.map(prepareAnnualPhenology);
      var modis = annual.mean();
      var greenup = modis.select('greenup').rename('Pheno_Greenup');
      var peak = modis.select('peak').rename('Pheno_Peak');
      var senescence = modis.select('senescence');
      var gsl = modis.select('gsl');
      var greenupDuration = modis.select('greenupDuration');
      var senescenceDuration = modis.select('senescenceDuration');
      var phenoAsymmetry = modis.select('asymmetry');
      var amp = modis.select('amp');

      // -----------------------------------------------------------------------
      // 4.4 Stack for neighborhood stats (use DOY for spatial variability)
      // -----------------------------------------------------------------------
      var phenoStack = ee.Image.cat([greenup, peak, gsl, amp])
        .rename(['greenup', 'peak', 'gsl', 'amp'])
        .reproject(ee.Image(modisCol.first()).select('Greenup_1').projection());
      
      // 1km scale
      var kernel_1km = ee.Kernel.circle({radius: CONFIG.modis.scales.small, units: 'meters'});
      var stats_1km = phenoStack.reduceNeighborhood({
        reducer: ee.Reducer.stdDev().combine(ee.Reducer.minMax(), null, true),
        kernel: kernel_1km,
        skipMasked: true
      });
      
      var greenup_sd_1km = stats_1km.select('greenup_stdDev').rename('Pheno_Greenup_SD_1km');
      var gsl_sd_1km = stats_1km.select('gsl_stdDev').rename('Pheno_GSL_SD_1km');
      var peak_sd_1km = stats_1km.select('peak_stdDev').rename('Pheno_Peak_SD_1km');
      var amp_sd_1km = stats_1km.select('amp_stdDev').rename('Pheno_Amp_SD_1km');
      var greenup_range_1km = stats_1km.select('greenup_max')
        .subtract(stats_1km.select('greenup_min')).rename('Pheno_Greenup_Range_1km');
      
      // 2km scale
      var kernel_2km = ee.Kernel.circle({radius: CONFIG.modis.scales.large, units: 'meters'});
      var stats_2km = phenoStack.select(['greenup', 'gsl']).reduceNeighborhood({
        reducer: ee.Reducer.stdDev(),
        kernel: kernel_2km,
        skipMasked: true
      });
      var greenup_sd_2km = stats_2km.select('greenup_stdDev').rename('Pheno_Greenup_SD_2km');
      var gsl_sd_2km = stats_2km.select('gsl_stdDev').rename('Pheno_GSL_SD_2km');
      
      // -----------------------------------------------------------------------
      // 4.5 Diversity indices
      // -----------------------------------------------------------------------
      var phenoDiversity_1km = greenup_sd_1km.add(gsl_sd_1km).add(peak_sd_1km)
        .divide(3).rename('Pheno_Diversity_1km');
      var phenoDiversity_2km = greenup_sd_2km.add(gsl_sd_2km)
        .divide(2).rename('Pheno_Diversity_2km');
      
      // -----------------------------------------------------------------------
      // 4.6 Stack all phenology features
      // -----------------------------------------------------------------------
      return ee.Image.cat([
        greenup,   // DOY
        peak,      // DOY
        senescence.rename('Pheno_Senescence'),  // DOY
        amp.rename('Pheno_Amp'),
        gsl.rename('Pheno_GSL'),                // Duration in days
        greenupDuration.rename('Pheno_GreenupDur'),      // Duration in days
        senescenceDuration.rename('Pheno_SenescenceDur'), // Duration in days
        phenoAsymmetry.rename('Pheno_Asymmetry'),         // Difference in days
        greenup_sd_1km, gsl_sd_1km, peak_sd_1km, amp_sd_1km, greenup_range_1km,
        greenup_sd_2km, gsl_sd_2km,
        phenoDiversity_1km, phenoDiversity_2km
      ]).toFloat();
    })(),
    ee.Image.constant(ee.List.repeat(0, outputBandNames.length))
      .rename(outputBandNames).toFloat().updateMask(ee.Image.constant(0))
  );
}

// ===================================================================================
// 6. MODIS Productivity Function (18 features)
// ===================================================================================
function getModisProductivity(year) {
  var outputBandNames = [
    // fPAR (6)
    'Prod_fPAR_Mean', 'Prod_fPAR_Max', 'Prod_fPAR_CV',
    'Prod_fPAR_SD_1km', 'Prod_fPAR_SD_2km', 'Prod_fPAR_Range_1km',
    // LAI (6)
    'Prod_LAI_Mean', 'Prod_LAI_Max', 'Prod_LAI_CV',
    'Prod_LAI_SD_1km', 'Prod_LAI_SD_2km', 'Prod_LAI_Range_1km',
    // GPP (6)
    'Prod_GPP_Mean', 'Prod_GPP_Sum', 'Prod_GPP_CV',
    'Prod_GPP_SD_1km', 'Prod_GPP_SD_2km', 'Prod_GPP_Range_1km'
  ];
  
  if (!CONFIG.productivity.enabled) {
    return ee.Image.constant(ee.List.repeat(0, outputBandNames.length))
      .rename(outputBandNames).toFloat().updateMask(ee.Image.constant(0));
  }
  
  // Three-year window: target year +/- 1 year.
  var yearNum = ee.Number(year);
  var startYear = yearNum.subtract(CONFIG.yearBufferForComposite);
  var endYear = yearNum.add(CONFIG.yearBufferForComposite);
  
  var startMonthStr = ee.Number(CONFIG.productivity.startMonth).format('%02d');
  var endMonthStr = ee.Number(CONFIG.productivity.endMonth).format('%02d');
  
  var start = ee.Date.parse('YYYY-MM-dd', 
    startYear.format('%d').cat('-').cat(startMonthStr).cat('-01'));
  var end = ee.Date.parse('YYYY-MM-dd', 
    endYear.format('%d').cat('-').cat(endMonthStr).cat('-30')).advance(1, 'day');
  
  // -------------------------------------------------------------------------
  // 5.1 fPAR and LAI from MOD15A2H (500m, 8-day) - 3 years of growing seasons
  // -------------------------------------------------------------------------
  var fpar_lai_col = ee.ImageCollection('MODIS/061/MOD15A2H')
    .filterDate(start, end)
    .map(maskMOD15A2H)
    // Growing-season months only.
    .filter(ee.Filter.calendarRange(
      CONFIG.productivity.startMonth, 
      CONFIG.productivity.endMonth, 
      'month'))
    .map(function(img) {
      var fpar = img.select('Fpar_500m').multiply(0.01);
      var lai = img.select('Lai_500m').multiply(0.1);
      
      fpar = fpar.updateMask(fpar.gte(0).and(fpar.lte(CONFIG.productivity.fparMax)));
      lai = lai.updateMask(lai.gte(0).and(lai.lte(CONFIG.productivity.laiMax)));
      
      return fpar.addBands(lai)
                 .rename(['fPAR', 'LAI'])
                 .copyProperties(img, ['system:time_start']);
    });
  
  var fpar_lai_size = fpar_lai_col.size();
  
  // -------------------------------------------------------------------------
  // 5.2 GPP from MOD17A2H (500m, 8-day) - 3 years of growing seasons
  // -------------------------------------------------------------------------
  var gpp_col = ee.ImageCollection('MODIS/061/MOD17A2H')
    .filterDate(start, end)
    .map(maskMOD17A2H)
    .filter(ee.Filter.calendarRange(
      CONFIG.productivity.startMonth, 
      CONFIG.productivity.endMonth, 
      'month'))
    .map(function(img) {
      var gpp = img.select('Gpp').multiply(0.0001);
      gpp = gpp.updateMask(gpp.gte(0).and(gpp.lte(CONFIG.productivity.gppMax)));
      
      return gpp.rename('GPP')
                .copyProperties(img, ['system:time_start']);
    });
  
  var gpp_size = gpp_col.size();
  
  // -------------------------------------------------------------------------
  // 5.3 Compute temporal statistics over 3-year window
  // -------------------------------------------------------------------------
  return ee.Algorithms.If(
    fpar_lai_size.gt(0).and(gpp_size.gt(0)),
    (function() {
      // fPAR statistics (3-year mean)
      var fpar_mean = fpar_lai_col.select('fPAR').mean().rename('Prod_fPAR_Mean');
      var fpar_max = fpar_lai_col.select('fPAR').max().rename('Prod_fPAR_Max');
      var fpar_stdDev = fpar_lai_col.select('fPAR').reduce(ee.Reducer.stdDev());
      var fpar_cv = fpar_stdDev.divide(fpar_mean.add(0.001)).rename('Prod_fPAR_CV');
      
      // LAI statistics (3-year mean)
      var lai_mean = fpar_lai_col.select('LAI').mean().rename('Prod_LAI_Mean');
      var lai_max = fpar_lai_col.select('LAI').max().rename('Prod_LAI_Max');
      var lai_stdDev = fpar_lai_col.select('LAI').reduce(ee.Reducer.stdDev());
      var lai_cv = lai_stdDev.divide(lai_mean.add(0.001)).rename('Prod_LAI_CV');
      
      // GPP statistics (3-year mean)
      var gpp_mean = gpp_col.mean().rename('Prod_GPP_Mean');
      // GPP_Sum: average annual growing season GPP (mean of yearly sums)
      var gpp_sum = gpp_col.sum().divide(CONFIG.yearBufferForComposite * 2 + 1)
        .rename('Prod_GPP_Sum');
      var gpp_stdDev = gpp_col.reduce(ee.Reducer.stdDev());
      var gpp_cv = gpp_stdDev.divide(gpp_mean.add(0.0001)).rename('Prod_GPP_CV');
      
      // -----------------------------------------------------------------------
      // 5.4 Compute spatial heterogeneity at 1km and 2km scales
      // -----------------------------------------------------------------------
      var prodStack = ee.Image.cat([fpar_mean, lai_mean, gpp_mean])
        .rename(['fpar', 'lai', 'gpp'])
        .reproject(ee.Image(fpar_lai_col.first()).select('fPAR').projection());
      
      // 1km scale
      var kernel_1km = ee.Kernel.circle({radius: CONFIG.modis.scales.small, units: 'meters'});
      var stats_1km = prodStack.reduceNeighborhood({
        reducer: ee.Reducer.stdDev().combine(ee.Reducer.minMax(), null, true),
        kernel: kernel_1km,
        skipMasked: true
      });
      
      var fpar_sd_1km = stats_1km.select('fpar_stdDev').rename('Prod_fPAR_SD_1km');
      var lai_sd_1km = stats_1km.select('lai_stdDev').rename('Prod_LAI_SD_1km');
      var gpp_sd_1km = stats_1km.select('gpp_stdDev').rename('Prod_GPP_SD_1km');
      
      var fpar_range_1km = stats_1km.select('fpar_max')
        .subtract(stats_1km.select('fpar_min')).rename('Prod_fPAR_Range_1km');
      var lai_range_1km = stats_1km.select('lai_max')
        .subtract(stats_1km.select('lai_min')).rename('Prod_LAI_Range_1km');
      var gpp_range_1km = stats_1km.select('gpp_max')
        .subtract(stats_1km.select('gpp_min')).rename('Prod_GPP_Range_1km');
      
      // 2km scale
      var kernel_2km = ee.Kernel.circle({radius: CONFIG.modis.scales.large, units: 'meters'});
      var stats_2km = prodStack.reduceNeighborhood({
        reducer: ee.Reducer.stdDev(),
        kernel: kernel_2km,
        skipMasked: true
      });
      
      var fpar_sd_2km = stats_2km.select('fpar_stdDev').rename('Prod_fPAR_SD_2km');
      var lai_sd_2km = stats_2km.select('lai_stdDev').rename('Prod_LAI_SD_2km');
      var gpp_sd_2km = stats_2km.select('gpp_stdDev').rename('Prod_GPP_SD_2km');
      
      // -----------------------------------------------------------------------
      // 5.5 Stack all productivity features
      // -----------------------------------------------------------------------
      return ee.Image.cat([
        // fPAR (6)
        fpar_mean, fpar_max, fpar_cv,
        fpar_sd_1km, fpar_sd_2km, fpar_range_1km,
        // LAI (6)
        lai_mean, lai_max, lai_cv,
        lai_sd_1km, lai_sd_2km, lai_range_1km,
        // GPP (6)
        gpp_mean, gpp_sum, gpp_cv,
        gpp_sd_1km, gpp_sd_2km, gpp_range_1km
      ]).toFloat();
    })(),
    ee.Image.constant(ee.List.repeat(0, outputBandNames.length))
      .rename(outputBandNames).toFloat().updateMask(ee.Image.constant(0))
  );
}

// ===================================================================================
// 7. Process One Year
// ===================================================================================
function processOneYear(targetYear) {
  var plotYear = ee.Number(targetYear);
  var plots_subset = plots.filter(ee.Filter.eq('year', plotYear));
  var offsets = CONFIG.augmentation.enabled ? CONFIG.augmentation.offsets : [0];
  
  var featuresForYear = offsets.map(function(offset) {
    var imageYear = plotYear.add(offset); 
    
    // Get phenology (17 features) - 3-year window with DOY
    var modisPheno = ee.Image(getModisPhenology(imageYear));
    
    // Get productivity (18 features) - 3-year window
    var modisProd = ee.Image(getModisProductivity(imageYear));
    
    // Combine all MODIS features (35 total)
    var modisAll = modisPheno.addBands(modisProd);
    
    return modisAll.reduceRegions({
      collection: plots_subset, 
      reducer: ee.Reducer.first(), 
      scale: CONFIG.export.scale, 
      crs: albers_wkt, 
      tileScale: CONFIG.tileScale
    }).map(function(f) {
      return f.set('plot_year', plotYear)
              .set('image_year', imageYear)
              .set('aug_offset', offset)
              .setGeometry(null);
    });
  });
  
  var training_features = ee.FeatureCollection(featuresForYear).flatten();
  training_features = training_features.filter(ee.Filter.notNull(['Pheno_Greenup', 'Prod_fPAR_Mean', 'Prod_LAI_Mean', 'Prod_GPP_Mean']));
  
  return training_features;
}

// ===================================================================================
// 8. Batch Export
// ===================================================================================
var allFeatures = ALL_YEARS.map(function(year) {
  return processOneYear(year);
});

var training_features = ee.FeatureCollection(allFeatures).flatten();

Export.table.toDrive({
  collection: training_features, 
  description: TEST_MODE ? 'MODIS_plot_predictors_TEST_2010' : 'MODIS_plot_predictors_AllYears', 
  folder: CONFIG.export.folder, 
  fileNamePrefix: TEST_MODE ? 'modis_plot_predictors_test' : 'modis_plot_predictors_all_years',
  fileFormat: 'CSV'
});

// ===================================================================================
// 9. Diagnostics
// ===================================================================================
print('============================================================');
print('MODIS productivity and phenology extraction');
print('============================================================');
print('Test mode: ' + (TEST_MODE ? 'ON (2010 only)' : 'OFF'));
print('Temporal window: target year +/- ' + CONFIG.yearBufferForComposite + ' year(s)');
print('Growing season: May-September');
print('MODIS sampling: nearest-neighbour/default assignment to the 30 m analysis grid');
print('Augmentation: ' + (CONFIG.augmentation.enabled ? CONFIG.augmentation.offsets.join(',') : 'OFF (aug_offset = 0 only)'));
print('Export folder: ' + CONFIG.export.folder);
print('Export prefix: ' + (TEST_MODE ? 'modis_plot_predictors_test' : 'modis_plot_predictors_all_years'));
print('Join keys for Landsat/MODIS tables: plot_id + year + aug_offset; plot_year and image_year are retained as audit fields');

var testYear = ALL_YEARS[0];
var testPlots = plots.filter(ee.Filter.eq('year', testYear));
print('Diagnostic plot count for ' + testYear + ':', testPlots.size());
print('Phenology bands:', ee.Image(getModisPhenology(testYear)).bandNames());
print('Productivity bands:', ee.Image(getModisProductivity(testYear)).bandNames());
