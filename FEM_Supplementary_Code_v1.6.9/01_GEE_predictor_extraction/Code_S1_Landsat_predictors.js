/**
 * Code S1. Landsat predictor extraction.
 *
 * Purpose
 *   Extract Landsat Collection 2 Level-2 surface-reflectance predictors,
 *   Landsat-derived spatial/temporal metrics, and static topographic metrics for model training.
 *
 * Method alignment
 *   - Target-year +/- 1 three-year temporal window.
 *   - Landsat 5 TM, Landsat 7 ETM+, and Landsat 8 OLI.
 *   - Landsat 8 OLI reflectance harmonized to the TM/ETM+ scale before pooling.
 *   - Spring = May-June, summer = July-August, autumn = September-October.
 *   - Spring/summer composites use a 90th-percentile-NDVI medoid.
 *   - Autumn composites use a median-NDVI medoid.
 *   - Extraction uses the Northeast China Albers equal-area projection at 30 m.
 */

// ===================================================================================
// 1. Global Configuration
// ===================================================================================
var CONFIG = {
  // Earth Engine feature collection containing plot coordinates and years.
  assetId: '', 
  
  // Keep aug_offset=0 downstream for backward-compatible merge keys.
  augmentation: {
    enabled: false,
    offsets: [0]
  },

  // Three-year window: target year +/- 1 year.
  yearBufferForComposite: 1,
  
  // Defined strategies for each season
  seasons: {
    spring: { monthStart: 5, monthEnd: 6, suffix: '_Spr', strategy: 'P90' },
    summer: { monthStart: 7, monthEnd: 8, suffix: '_Sum', strategy: 'P90' },
    autumn: { monthStart: 9, monthEnd: 10, suffix: '_Aut', strategy: 'Median' }
  },
  
  texture: {
    enabled: true,
    glcmSize: 3,
    metrics: ['contrast', 'ent', 'asm', 'corr', 'var', 'diss'],
    targetBands: ['NIR', 'SWIR1', 'EVI', 'NDVI']
  },
  
  entropy: {
    enabled: true,
    kernelRadius: 3,
    targetBands: ['NIR', 'SWIR1', 'NDVI', 'EVI', 'NDMI']
  },
  
  spectralMorphology: {
    enabled: true,
    targetBands: ['NDVI', 'EVI', 'NDMI', 'NIR', 'SWIR1']
  },
stms: {
    enabled: true,
    targetBands: ['NDVI', 'EVI', 'SAVI', 'NBR', 'NDMI', 'VARI', 'NDSI', 'GRVI', 'NIR', 'SWIR1', 'SWIR2']
  },
  
  spectralAngle: { 
    enabled: true, 
    bands: ['B1', 'B2', 'B3', 'NIR', 'SWIR1', 'SWIR2'] 
  },
  
  multiScale: { windowSizes: [1, 2, 3] },
  
  largeScale: {
    radiiMeters: [450, 900, 1800],
    targetBandsCV: ['NDVI', 'EVI', 'NDMI', 'Elev', 'Slope', 'TPI_1km']
  },
  
  smallScaleBands: ['B1', 'B2', 'B3', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'EVI', 'SAVI', 'NBR', 'NDMI', 'VARI', 'NDSI', 'GRVI'],

  seasonalBands: ['B1', 'B2', 'B3', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'EVI', 'SAVI', 'NBR', 'NDMI', 'VARI', 'NDSI', 'GRVI'],

  diffBands: ['NDVI', 'EVI', 'SAVI', 'NBR', 'NDMI', 'VARI', 'NDSI', 'GRVI', 'NIR', 'SWIR1'],
  
  export: { folder: 'plot_predictor_exports', scale: 30 },
  
  tileScale: 16, 
  epsilon: 1e-6
};

// Optional test mode
var TEST_MODE = false;
var ALL_YEARS = TEST_MODE ? [2010] : [2008, 2009, 2010, 2011, 2012, 2014, 2015, 2016, 2017];

// ===================================================================================
// 2. Projection & Ancillary Data
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

var TOPO_BASE = (function() {
  var nasadem = ee.Image('NASA/NASADEM_HGT/001').select('elevation').toFloat().rename('Elev');
  var slope = ee.Terrain.slope(nasadem).rename('Slope');
  var aspect = ee.Terrain.aspect(nasadem);
  var northness = aspect.multiply(Math.PI / 180).cos().rename('Northness');
  var eastness = aspect.multiply(Math.PI / 180).sin().rename('Eastness');
  var tpi_1km = nasadem.subtract(nasadem.focalMean(1000, 'circle', 'meters')).rename('TPI_1km');
  var tri = nasadem.reduceNeighborhood({ 
    reducer: ee.Reducer.stdDev(), 
    kernel: ee.Kernel.square(1) 
  }).rename('TRI');
  return nasadem.addBands([slope, northness, eastness, tpi_1km, tri]);
})();

// ===================================================================================
// 3. Image Pre-processing Functions
// ===================================================================================
function maskClouds(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));
  return image.updateMask(mask);
}

function applyScaleFactors(image) {
  var sr = image.select('SR_B.');

  // Landsat Collection 2 SR valid range: DN 7273–43636 maps approximately to 0–1 reflectance.
  var valid = sr.reduce(ee.Reducer.min()).gte(7273)
    .and(sr.reduce(ee.Reducer.max()).lte(43636));

  var scaled = sr.multiply(0.0000275).add(-0.2).updateMask(valid);
  return image.addBands(scaled, null, true);
}

function harmonizeOLI(image) {
  var bands = image.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']);

  // OLI -> ETM+ OLS surface-reflectance harmonization from Roy et al. (2016), Table 2.
  // Band order: blue, green, red, NIR, SWIR1, SWIR2.
  var slopes = ee.Image.constant([0.8850, 0.9317, 0.9372, 0.8339, 0.8639, 0.9165]);
  var intercepts = ee.Image.constant([0.0183, 0.0123, 0.0123, 0.0448, 0.0306, 0.0116]);

  return image.addBands(bands.multiply(slopes).add(intercepts), null, true);
}

function renameL8(image) { 
  return image.select(
    ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'], 
    ['B1','B2','B3','NIR','SWIR1','SWIR2']
  ); 
}

function renameL57(image) { 
  return image.select(
    ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'], 
    ['B1','B2','B3','NIR','SWIR1','SWIR2']
  ); 
}

function addAllIndices(image) {
  var n = image.select('NIR');
  var r = image.select('B3');  // Red
  var g = image.select('B2');  // Green
  var b = image.select('B1');  // Blue
  var s1 = image.select('SWIR1');
  var s2 = image.select('SWIR2');
  var eps = CONFIG.epsilon;
  
  var ndvi = n.subtract(r).divide(n.add(r).add(eps)).rename('NDVI');
  var evi = image.expression('2.5*((N-R)/(N+6*R-7.5*B+1))',{'N':n,'R':r,'B':b}).rename('EVI');
  var savi = image.expression('1.5*((N-R)/(N+R+0.5))',{'N':n,'R':r}).rename('SAVI');
  var nbr = n.subtract(s2).divide(n.add(s2).add(eps)).rename('NBR');
  var ndmi = n.subtract(s1).divide(n.add(s1).add(eps)).rename('NDMI');
  var ndsi = g.subtract(s1).divide(g.add(s1).add(eps)).rename('NDSI');
  
  var grvi = g.subtract(r).divide(g.add(r).add(eps)).rename('GRVI');
  var variDenom = g.add(r).subtract(b);
  var vari = g.subtract(r).divide(variDenom.where(variDenom.abs().lt(eps), eps)).rename('VARI');
  
  return image.addBands([ndvi, evi, savi, nbr, ndmi, ndsi, grvi, vari]).toFloat();
}

function preprocessL57(image) {
  return addAllIndices(renameL57(maskClouds(applyScaleFactors(image))));
}

function preprocessL89(image) {
  return addAllIndices(renameL8(harmonizeOLI(maskClouds(applyScaleFactors(image)))));
}


// ===================================================================================
// 5. Entropy Calculation
// ===================================================================================
function computeEntropy(image, suffix) {
  if (!CONFIG.entropy.enabled) {
    return ee.Image([]);
  }
  
  var kernel = ee.Kernel.square({radius: CONFIG.entropy.kernelRadius});
  var targetBands = CONFIG.entropy.targetBands;
  
  var bandList = targetBands.map(function(band) {
    var bandImg = image.select(band);
    var isIndex = ['NDVI', 'EVI', 'SAVI', 'NBR', 'NDMI'].indexOf(band) > -1;
    
    var scaled = isIndex 
      ? bandImg.clamp(-1, 1).add(1).multiply(127.5).toByte()
      : bandImg.clamp(0, 1).multiply(255).toByte();
    
    return scaled.entropy(kernel).rename(band + '_Ent' + suffix);
  });
  
  return ee.Image.cat(bandList);
}

// ===================================================================================
// 6. Spectral Morphology
// ===================================================================================
function computeSpectralMorphology(image, suffix) {
  if (!CONFIG.spectralMorphology.enabled) {
    return ee.Image([]);
  }
  
  var targetBands = CONFIG.spectralMorphology.targetBands;
  var selectedImage = image.select(targetBands);
  
  var dilation = selectedImage.spectralDilation('sam');
  var erosion = selectedImage.spectralErosion('sam');
  var gradient = dilation.subtract(erosion);
  
  var dilationRenamed = dilation.rename(targetBands.map(function(b) { return b + '_SDil' + suffix; }));
  var erosionRenamed = erosion.rename(targetBands.map(function(b) { return b + '_SEro' + suffix; }));
  var gradientRenamed = gradient.rename(targetBands.map(function(b) { return b + '_SGrad' + suffix; }));
  
  return dilationRenamed.addBands(erosionRenamed).addBands(gradientRenamed);
}

// ===================================================================================
// 7. STMs - Spectral Temporal Metrics
// ===================================================================================
function calculateSTMs(collection, suffix) {
  if (!CONFIG.stms.enabled) {
    return ee.Image([]);
  }
  
  var targetBands = CONFIG.stms.targetBands;
  
  var combinedReducer = ee.Reducer.mean()
    .combine(ee.Reducer.max(), null, true)
    .combine(ee.Reducer.min(), null, true)
    .combine(ee.Reducer.stdDev(), null, true)
    .combine(ee.Reducer.percentile([10, 50, 90]), null, true);
  
  var stats = collection.select(targetBands).reduce(combinedReducer);
  
  var renamedBands = [];
  targetBands.forEach(function(band) {
    renamedBands.push(stats.select(band + '_mean').rename(band + '_Mean' + suffix));
    renamedBands.push(stats.select(band + '_max').rename(band + '_Max' + suffix));
    renamedBands.push(stats.select(band + '_min').rename(band + '_Min' + suffix));
    renamedBands.push(stats.select(band + '_stdDev').rename(band + '_Std' + suffix));
    renamedBands.push(stats.select(band + '_p10').rename(band + '_P10' + suffix));
    renamedBands.push(stats.select(band + '_p50').rename(band + '_Med' + suffix));
    renamedBands.push(stats.select(band + '_p90').rename(band + '_P90' + suffix));
  });
  
  var derivedBands = targetBands.map(function(band) {
    return stats.select(band + '_max').subtract(stats.select(band + '_min'))
      .rename(band + '_Range' + suffix);
  });
  
  return ee.Image.cat(renamedBands).addBands(ee.Image.cat(derivedBands));
}

// ===================================================================================
// 8. Seasonal Collection Builder
// ===================================================================================
function getSeasonalCollection(imageYear, geom, seasonConfig) {
  var year = ee.Number(imageYear);
  var startYear = year.subtract(CONFIG.yearBufferForComposite);
  var endYear = year.add(CONFIG.yearBufferForComposite);
  
  var start_date = ee.Date.fromYMD(startYear, 1, 1);
  var end_date = ee.Date.fromYMD(endYear, 12, 31);
  var monthFilter = ee.Filter.calendarRange(seasonConfig.monthStart, seasonConfig.monthEnd, 'month');
  
  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(geom).filterDate(start_date, end_date).filter(monthFilter).map(preprocessL57);
  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(geom).filterDate(start_date, end_date).filter(monthFilter).map(preprocessL57);
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(geom).filterDate(start_date, end_date).filter(monthFilter).map(preprocessL89);
  
  // Study period and method use Landsat 5, 7, and 8 only.
  return l5.merge(l7).merge(l8);
}

// ===================================================================================
// 9. Texture Computation
// ===================================================================================
function computeTexture(image, kernel, suffix) {
  var bands = CONFIG.texture.targetBands;
  var metrics = CONFIG.texture.metrics;
  
  var textureInputs = ee.Image.cat(bands.map(function(bandName) {
    var input = image.select(bandName);
    var isIndex = ['EVI','NDVI','SAVI','NBR','NDMI'].indexOf(bandName) > -1;
    return isIndex 
      ? input.clamp(-1, 1).add(1).multiply(10000).toInt16()
      : input.clamp(0, 2.0).multiply(10000).toInt16();
  })).rename(bands);
  
  var glcm = textureInputs.glcmTexture({size: CONFIG.texture.glcmSize});
  
  var selectedBands = [];
  bands.forEach(function(b) {
    metrics.forEach(function(m) {
      selectedBands.push(b + '_' + m);
    });
  });
  
  var smoothed = glcm.select(selectedBands)
    .reduceNeighborhood({ reducer: ee.Reducer.mean(), kernel: kernel });
  
  var newNames = smoothed.bandNames().map(function(n) { 
    return ee.String(n).cat(suffix); 
  });
  
  return smoothed.rename(newNames);
}

// ===================================================================================
// 10. Spectral Indices (SG, SD)
// ===================================================================================
function computeSpectralIndices(image, kernel, suffix) {
  var bands = CONFIG.spectralAngle.bands;
  
  var combinedReducer = ee.Reducer.stdDev()
    .combine(ee.Reducer.minMax(), null, true);
  
  var stats = image.select(bands).reduceNeighborhood({
    reducer: combinedReducer, 
    kernel: kernel
  });
  
  var stdDevBands = bands.map(function(b) { return b + '_stdDev'; });
  var sg = stats.select(stdDevBands).pow(2).reduce(ee.Reducer.mean()).sqrt().rename('SG' + suffix);
  
  var minBands = bands.map(function(b) { return b + '_min'; });
  var maxBands = bands.map(function(b) { return b + '_max'; });
  var ranges = stats.select(maxBands).subtract(stats.select(minBands));
  var sd = ranges.pow(2).reduce(ee.Reducer.mean()).sqrt().rename('SD' + suffix);
  
  return sg.addBands(sd);
}

// ===================================================================================
// 11. Calculate Features at Scale
// ===================================================================================
function calculateFeaturesAtScale(image, radius, includeTextureAndEntropy, seasonSuffix) {
  var kernel = ee.Kernel.square(radius, 'pixels');
  var scaleSuffix = '_w' + (radius * 2 + 1);
  var suffix = scaleSuffix + seasonSuffix;
  var targetBands = CONFIG.smallScaleBands;
  
  var stats = image.select(targetBands).reduceNeighborhood({
    reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true),
    kernel: kernel
  });
  
  var statsRenamed = stats.rename(stats.bandNames().map(function(n) { 
    return ee.String(n).cat(suffix); 
  }));
  
  var cvList = targetBands.map(function(b) {
    var meanBand = b + '_mean' + suffix;
    var sdBand = b + '_stdDev' + suffix;
    return statsRenamed.select(sdBand)
      .divide(statsRenamed.select(meanBand).abs().add(CONFIG.epsilon))
      .rename(b + '_CV' + suffix);
  });
  var cvImage = ee.Image.cat(cvList);
  
  var result = statsRenamed.addBands(cvImage);
  
  if (includeTextureAndEntropy && radius <= 2) {
    if (CONFIG.texture.enabled) {
      result = result.addBands(computeTexture(image, kernel, suffix));
    }
    if (CONFIG.entropy.enabled) {
      result = result.addBands(computeEntropy(image, suffix));
    }
  }
  
  if (CONFIG.spectralAngle.enabled) {
    result = result.addBands(computeSpectralIndices(image, kernel, suffix));
  }
  
  return result;
}

// ===================================================================================
// 12. Large Scale Features
// ===================================================================================
function calculateLargeScaleFeatures(image, seasonSuffix) {
  var targetBandsCV = CONFIG.largeScale.targetBandsCV;
  var scalesMeters = CONFIG.largeScale.radiiMeters;
  var spectralBands = CONFIG.spectralAngle.bands;
  
  var scaleImages = scalesMeters.map(function(m) {
    var kernel = ee.Kernel.square(m, 'meters');
    var suffix = '_' + m + 'm' + seasonSuffix;
    
    var stats = image.select(targetBandsCV).reduceNeighborhood({
      reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true), 
      kernel: kernel
    });
    
    var cvList = targetBandsCV.map(function(b) {
      return stats.select(b + '_stdDev')
        .divide(stats.select(b + '_mean').abs().add(CONFIG.epsilon))
        .rename(b + '_CV' + suffix);
    });
    var resultBlock = ee.Image.cat(cvList);
    
    if (CONFIG.spectralAngle.enabled) {
      var spectralStats = image.select(spectralBands).reduceNeighborhood({
        reducer: ee.Reducer.stdDev(),
        kernel: kernel
      });
      var aSG = spectralStats.pow(2).reduce(ee.Reducer.mean()).sqrt().rename('aSG' + suffix);
      resultBlock = resultBlock.addBands(aSG);
    }
    
    return resultBlock;
  });
  
  return ee.ImageCollection.fromImages(scaleImages).toBands().regexpRename('^[0-9]+_', '');
}


// ===================================================================================
// 13. Seasonal Medoid Compositing
// ===================================================================================
function buildSeasonalComposite(col, seasonSuffix, strategy) {
  var bandsList = CONFIG.seasonalBands;
  var medoidImg;

  if (strategy === 'Median') {
    // Autumn: select the observation closest to seasonal median NDVI.
    var ndviTarget = col.select('NDVI').reduce(ee.Reducer.median()).unmask(-1);
    
    var colWithDist = col.map(function(img) {
      var dist = img.select('NDVI').subtract(ndviTarget).abs().rename('inv_dist');
      // Invert distance so qualityMosaic() selects the smallest distance.
      return img.addBands(dist.multiply(-1));
    });
    
    medoidImg = colWithDist.qualityMosaic('inv_dist').select(bandsList);

  } else {
    // Spring and summer: select the observation closest to the seasonal 90th-percentile NDVI.
    var ndviTarget = col.select('NDVI').reduce(ee.Reducer.percentile([90])).unmask(-1);
    
    var colWithDist = col.map(function(img) {
      var dist = img.select('NDVI').subtract(ndviTarget).abs().rename('inv_dist');
      // Invert distance so qualityMosaic() selects the smallest distance.
      return img.addBands(dist.multiply(-1));
    });
    
    medoidImg = colWithDist.qualityMosaic('inv_dist').select(bandsList);
  }
  // Freeze the 30-m analysis support before any pixel-neighborhood operation.
  medoidImg = medoidImg.reproject({crs: albers_wkt, scale: CONFIG.export.scale});
  // Seasonal median and temporal-standard-deviation features.
  var simpleMedian = col.select(bandsList).median();
  
  var sdComposite = col.select(bandsList).reduce(ee.Reducer.stdDev());
  
  var medoidRenamed = medoidImg.rename(bandsList.map(function(n) { return n + seasonSuffix; }));
  var medianRenamed = simpleMedian.rename(bandsList.map(function(n) { return n + '_Med' + seasonSuffix; }));
  var sdRenamed = sdComposite.rename(bandsList.map(function(n) { return n + '_TSD' + seasonSuffix; }));
  
  return {
    peak: medoidImg,
    peakRenamed: medoidRenamed,
    median: medianRenamed,
    tempSD: sdRenamed
  };
}

// ===================================================================================
// 14. Calculate Seasonal Differences
// ===================================================================================
function calculateSeasonalDifferences(summerPeak, springPeak, autumnPeak) {
  var diffBands = CONFIG.diffBands;
  
  var sumSprDiff = summerPeak.select(diffBands).subtract(springPeak.select(diffBands))
    .rename(diffBands.map(function(b) { return b + '_SumSprDiff'; }));
  
  var autSumDiff = autumnPeak.select(diffBands).subtract(summerPeak.select(diffBands))
    .rename(diffBands.map(function(b) { return b + '_AutSumDiff'; }));
  
  var minSprAut = springPeak.select(diffBands).min(autumnPeak.select(diffBands));
  var amplitude = summerPeak.select(diffBands).subtract(minSprAut)
    .rename(diffBands.map(function(b) { return b + '_SeasonAmp'; }));
  
  return sumSprDiff.addBands(autSumDiff).addBands(amplitude);
}

// ===================================================================================
// 15. Core Construction
// ===================================================================================
function buildMultiSeasonFeatureStack(imageYear, geom) {
  
  var springCol = getSeasonalCollection(imageYear, geom, CONFIG.seasons.spring);
  var summerCol = getSeasonalCollection(imageYear, geom, CONFIG.seasons.summer);
  var autumnCol = getSeasonalCollection(imageYear, geom, CONFIG.seasons.autumn);
  
  var hasData = springCol.size().gte(1)
    .and(summerCol.size().gte(1))
    .and(autumnCol.size().gte(1));
  
  return ee.Algorithms.If(
    hasData,
    (function() {
      var summerColFilled = summerCol;
      
      var spring = buildSeasonalComposite(springCol, CONFIG.seasons.spring.suffix, CONFIG.seasons.spring.strategy);
      var summer = buildSeasonalComposite(summerColFilled, CONFIG.seasons.summer.suffix, CONFIG.seasons.summer.strategy);
      var autumn = buildSeasonalComposite(autumnCol, CONFIG.seasons.autumn.suffix, CONFIG.seasons.autumn.strategy);
      
      // Small Scale Features
      var springSmallScaleList = CONFIG.multiScale.windowSizes.map(function(r) { 
        return calculateFeaturesAtScale(spring.peak, r, false, CONFIG.seasons.spring.suffix); 
      });
      var springSmallScale = ee.ImageCollection.fromImages(springSmallScaleList)
        .toBands().regexpRename('^[0-9]+_', '');
      
      var summerSmallScaleList = CONFIG.multiScale.windowSizes.map(function(r) { 
        return calculateFeaturesAtScale(summer.peak, r, true, CONFIG.seasons.summer.suffix); 
      });
      var summerSmallScale = ee.ImageCollection.fromImages(summerSmallScaleList)
        .toBands().regexpRename('^[0-9]+_', '');
      
      var autumnSmallScaleList = CONFIG.multiScale.windowSizes.map(function(r) { 
        return calculateFeaturesAtScale(autumn.peak, r, false, CONFIG.seasons.autumn.suffix); 
      });
      var autumnSmallScale = ee.ImageCollection.fromImages(autumnSmallScaleList)
        .toBands().regexpRename('^[0-9]+_', '');
      
      // Large Scale Features
      var combinedSummer = summer.peak.addBands(TOPO_BASE);
      var largeScaleFeatures = calculateLargeScaleFeatures(combinedSummer, CONFIG.seasons.summer.suffix);
      
      // Seasonal Differences
      var seasonalDiffs = calculateSeasonalDifferences(summer.peak, spring.peak, autumn.peak);
      
      // Spectral Morphology
      var summerMorphology = computeSpectralMorphology(summer.peak, CONFIG.seasons.summer.suffix);
      
      // STMs
      var fullSeasonCol = springCol.merge(summerColFilled).merge(autumnCol);
      var stmsFeatures = calculateSTMs(fullSeasonCol, '_GS');
      
      // Merge All
      var featureStack = ee.Image([])
        .addBands(TOPO_BASE)
        .addBands(spring.peakRenamed)
        .addBands(summer.peakRenamed)
        .addBands(autumn.peakRenamed)
        .addBands(spring.median)
        .addBands(summer.median)
        .addBands(autumn.median)
        .addBands(spring.tempSD)
        .addBands(summer.tempSD)
        .addBands(autumn.tempSD)
        .addBands(springSmallScale)
        .addBands(summerSmallScale)
        .addBands(autumnSmallScale)
        .addBands(largeScaleFeatures)
        .addBands(seasonalDiffs)
        .addBands(summerMorphology)
        .addBands(stmsFeatures);
      
      return featureStack.toFloat();
    })(),
    ee.Image([])
  );
}

// ===================================================================================
// 16. Single Year Processing
// ===================================================================================
function processOneYear(targetYear) {
  var plotYear = ee.Number(targetYear);
  var plots_subset = plots.filter(ee.Filter.eq('year', plotYear));
  var offsets = CONFIG.augmentation.enabled ? CONFIG.augmentation.offsets : [0];
  
  var featuresForYear = offsets.map(function(offset) {
    var imageYear = plotYear.add(offset); 
    var featureStack = ee.Image(buildMultiSeasonFeatureStack(imageYear, plots_subset));
    
    return ee.Algorithms.If(
      featureStack.bandNames().length().gt(0),
      featureStack.reduceRegions({
        collection: plots_subset, 
        reducer: ee.Reducer.first(), 
        scale: CONFIG.export.scale, 
        crs: albers_wkt, 
        tileScale: CONFIG.tileScale
      }).map(function(f) {
        return f.set('status', 'SUCCESS')
                .set('plot_year', plotYear)
                .set('image_year', imageYear)
                .set('aug_offset', offset);
      }),
      plots_subset.map(function(f) {
        return f.set('status', 'FAILED_NO_DATA')
                .set('plot_year', plotYear)
                .set('image_year', imageYear);
      })
    );
  });
  
  var training_features = ee.FeatureCollection(featuresForYear).flatten();
  training_features = training_features.filter(ee.Filter.eq('status', 'SUCCESS'));
  
  return training_features;
}

// ===================================================================================
// 17. Batch Export
// ===================================================================================
var allFeatures = ALL_YEARS.map(function(year) {
  return processOneYear(year);
});

var training_features = ee.FeatureCollection(allFeatures).flatten();

Export.table.toDrive({
  collection: training_features, 
  description: TEST_MODE ? 'Landsat_plot_predictors_TEST' : 'Landsat_plot_predictors_AllYears', 
  folder: CONFIG.export.folder, 
  fileNamePrefix: TEST_MODE ? 'landsat_plot_predictors_test' : 'landsat_plot_predictors_all_years',
  fileFormat: 'CSV'
});

// ===================================================================================
// 18. Diagnostics
// ===================================================================================
print('============================================================');
print('LANDSAT predictor extraction');
print('============================================================');
print('Test mode: ' + (TEST_MODE ? 'ON' : 'OFF'));
print('Temporal window: target year +/- ' + CONFIG.yearBufferForComposite + ' year(s)');
print('Seasons: spring May-Jun, summer Jul-Aug, autumn Sep-Oct');
print('Composites: spring/summer P90-NDVI medoid; autumn median-NDVI medoid');
print('Augmentation: ' + (CONFIG.augmentation.enabled ? CONFIG.augmentation.offsets.join(',') : 'OFF (aug_offset = 0 only)'));
print('Export folder: ' + CONFIG.export.folder);
print('Export prefix: ' + (TEST_MODE ? 'landsat_plot_predictors_test' : 'landsat_plot_predictors_all_years'));
print('Spectral indices retained in the Landsat export: NDVI, EVI, SAVI, NBR, NDMI, VARI, NDSI, GRVI');
print('============================================================');