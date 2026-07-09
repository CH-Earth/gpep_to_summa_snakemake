// MODIS snow metrics by subbasin (CSV export)
// Matches compute_water_year_snow_metrics.py defaults where noted:
//   - SCF = snow days / valid days (snow counted only on valid observations)
//   - SDD = 5-day snow run followed by 5-day snow-free run (latest transition)
//
// Processes water years 2008–2019 using one reduceRegions call per export.

//////////////////////////////////////////////
////////   USER INPUT REQUIRED  //////////////
//////////////////////////////////////////////

var startWyr = 2008;
var endWyr = 2019;

var subbasins = ee.FeatureCollection('users/davecasson/chena_tdx');
var subbasinIdField = 'fid';

var waterYearStartMonth = 10;
var waterYearStartDay = 1;

var NDSI_threshold = 15;
var runDays = 5;

// Minimum clear MODIS observations at a pixel to report SDD.
// Python uses 360 on hourly model forcing; MODIS pixels rarely reach that.
// Set to 0 to disable this pixel filter.
var minValidWyDays = 100;

// Mask SCF (not SDD) where annual snow frequency is below this value.
// Set to 0 to disable.
var lowScfMask = 0.07;

var exportScale = 500;
var projection = 'EPSG:4326';
var exportFolder = 'GEE_snow_metrics_chena_tdx_csv';

var modisSnowCollection = 'MODIS/061/MOD10A1';
var region = ee.Geometry.Rectangle([-148, 64, -143, 66]);

//////////////////////////////////////////////
////////    END USER INPUT      //////////////
//////////////////////////////////////////////

print('Script started');
print('Processing water years', startWyr + ' to ' + endWyr);


/**
 * Daily 0/1 flags from MOD10A1.
 *
 * Valid = NDSI_Snow_Cover in [0, 100] (GEE masks class codes >100 in this band).
 * Do NOT use class != 250 alone; missing/night/no-decision classes are != 250
 * but are not usable observations and were breaking SDD persistence windows.
 */
var dailyFlagsFromModis = function(img, wyStart) {
  var snowBand = img.select('NDSI_Snow_Cover');

  var valid_v = snowBand
    .gte(0)
    .and(snowBand.lte(100))
    .toFloat()
    .unmask(0)
    .rename('valid_v');

  var snow_v = snowBand
    .gte(NDSI_threshold)
    .and(snowBand.lte(100))
    .toFloat()
    .unmask(0)
    .rename('snow_v');

  var nosnow_v = valid_v.multiply(snow_v.eq(0)).rename('nosnow_v');

  var doy = img
    .date()
    .difference(wyStart, 'day')
    .add(1)
    .rename('doy')
    .toFloat();

  return img.addBands([valid_v, snow_v, nosnow_v, doy])
    .select(['valid_v', 'snow_v', 'nosnow_v', 'doy'])
    .clip(region);
};


/**
 * SDD via runDays persistence rule using filterDate windows (no toList/list.get).
 *
 * For each day D (first snow-free day of the melt-out run):
 *   prior  runDays days: [D-runDays, D)     all valid and snow
 *   after  runDays days: [D, D+runDays)     all valid and snow-free
 *
 * Keep the latest qualifying transition (max day-of-water-year).
 */
var sddFromPersistence = function(dailyCol, runDays, minValidWyDays) {
  dailyCol = dailyCol.sort('system:time_start');
  var run = ee.Number(runDays);

  var transitions = dailyCol.map(function(img) {
    var date = img.date();

    var prior = dailyCol.filterDate(
      date.advance(run.multiply(-1), 'day'),
      date
    );
    var after = dailyCol.filterDate(
      date,
      date.advance(run, 'day')
    );

    var sumPriorSnow = prior.select('snow_v').sum();
    var sumPriorValid = prior.select('valid_v').sum();
    var sumAfterNoSnow = after.select('nosnow_v').sum();
    var sumAfterValid = after.select('valid_v').sum();

    var isTransition = sumPriorSnow.eq(run)
      .and(sumPriorValid.eq(run))
      .and(sumAfterNoSnow.eq(run))
      .and(sumAfterValid.eq(run));

    return img.select('doy').updateMask(isTransition).rename('candidate_doy');
  });

  var validCount = dailyCol.select('valid_v').sum();
  var sdd = transitions.select('candidate_doy').max();

  sdd = sdd.updateMask(sdd.gt(0));

  if (minValidWyDays > 0) {
    sdd = sdd.updateMask(validCount.gte(minValidWyDays));
  }

  return sdd.rename('sdd');
};


var getMetricImageForYear = function(wyr) {
  wyr = ee.Number(wyr);
  var wyrString = wyr.format('%d');

  var wyStart = ee.Date.fromYMD(
    wyr.subtract(1),
    waterYearStartMonth,
    waterYearStartDay
  );
  var wyEnd = ee.Date.fromYMD(
    wyr,
    waterYearStartMonth,
    waterYearStartDay
  );

  var modis = ee.ImageCollection(modisSnowCollection)
    .filterDate(wyStart, wyEnd)
    .filterBounds(region)
    .map(function(img) {
      return dailyFlagsFromModis(img, wyStart);
    });

  var snowCount = modis.select('snow_v').sum();
  var validCount = modis.select('valid_v').sum();

  var SCF = snowCount
    .divide(validCount)
    .rename(ee.String('SCF_').cat(wyrString))
    .updateMask(validCount.gt(0));

  if (lowScfMask > 0) {
    SCF = SCF.updateMask(SCF.gte(lowScfMask));
  }

  var SDD = sddFromPersistence(modis, runDays, minValidWyDays)
    .rename(ee.String('SDD_').cat(wyrString));

  // SDD is independent of the low-SCF mask; only require some valid obs in the pixel.
  SDD = SDD.updateMask(validCount.gt(0));

  return SCF.addBands(SDD).toFloat();
};


var years = ee.List.sequence(startWyr, endWyr);
var emptyImage = ee.Image([]);

var metricsImage = ee.Image(
  years.iterate(function(y, imageSoFar) {
    imageSoFar = ee.Image(imageSoFar);
    return imageSoFar.addBands(getMetricImageForYear(y));
  }, emptyImage)
);

print('Metric image bands', metricsImage.bandNames());


// --- Optional spot-check (uncomment before export) ---
// var test2015 = getMetricImageForYear(2015);
// var testGeom = subbasins.first().geometry();
// print('2015 HRU means', test2015.reduceRegion({
//   reducer: ee.Reducer.mean(),
//   geometry: testGeom,
//   scale: exportScale,
//   maxPixels: 1e9
// }));
// Map.addLayer(test2015.select('SCF_2015'), {min: 0, max: 1}, 'SCF 2015');
// Map.addLayer(test2015.select('SDD_2015'), {min: 150, max: 200}, 'SDD 2015');


var zonalReducer = ee.Reducer.mean()
  .combine({reducer2: ee.Reducer.median(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.minMax(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.count(), sharedInputs: true});

var statsWideRaw = metricsImage.reduceRegions({
  collection: subbasins,
  reducer: zonalReducer,
  scale: exportScale,
  crs: projection,
  tileScale: 8
});

var statsWide = statsWideRaw.map(function(f) {
  var props = f.toDictionary();
  props = props.set('fid', f.get(subbasinIdField));
  return ee.Feature(null, props);
});

print('Sample output', statsWide.limit(5));

Export.table.toDrive({
  collection: statsWide,
  description: 'SnowMetrics_2008_2019_chena_tdx_subbasin_stats_wide',
  folder: exportFolder,
  fileNamePrefix: 'SnowMetrics_2008_2019_chena_tdx_subbasin_stats_wide',
  fileFormat: 'CSV'
});

print('Export task created for water years 2008–2019');
