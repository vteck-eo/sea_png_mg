/************************************************************
 🌿 Mangrove Boundary Extraction from Probability Surface
************************************************************/


/************************************************************
 LOAD EMBEDDING-BASED MANGROVE PROBABILITY (10 m)
*************************************************************/

var assetList = {
  2017: "image/2017"
};

function loadProb(year) {
  return ee.Image(assetList[year])
    .select("classification")
    .divide(100)                 // convert 0–100 → 0–1
    .rename("p" + year);
}

var p = loadProb(2017);

Map.addLayer(
  p,
  {min: 0, max: 1, palette: ['blue', 'yellow', 'red']},
  "Mangrove Probability 2017",
  false
);


/************************************************************
 DEFINE COUNTRY AOI (FAO GAUL 2015)
*************************************************************/

var fao = ee.FeatureCollection("FAO/GAUL/2015/level0");

var AOI = fao
  .filter(ee.Filter.eq('ADM0_NAME', 'Country'))
  .geometry()
  .buffer(5000);   // small buffer to avoid edge clipping

Map.addLayer(AOI, {}, "Boundary", false);

var aoi = AOI;

Map.addLayer(aoi, {color: 'cyan'}, 'AOI', false);


/************************************************************
 CREATE MANGROVE MASK (10 m Raster)
*************************************************************/

// Probability threshold (adjust if needed)
var threshold = 0.50;

var mangroveMask = p
  .gte(threshold)
  .selfMask()
  .clip(aoi);

Map.addLayer(
  mangroveMask,
  {palette: ['00FF00']},
  "Mangrove Mask (Thresholded)",
  false
);


/************************************************************
 RASTER BUFFER (Distance Transform Method)
*************************************************************/

// Buffer distance in meters
var BUFFER_M = 100;

// Distance transform (convert pixel distance → meters)
var dist = mangroveMask
  .fastDistanceTransform(30)
  .sqrt()
  .multiply(10);   // 10 m pixel size

var mangroveBuffered = dist
  .lte(BUFFER_M)
  .selfMask();

Map.addLayer(
  mangroveBuffered,
  {palette: ['FF0000']},
  "Mangrove Buffered (Raster)",
  false
);


/************************************************************
 APPLY MINIMUM MAPPING UNIT (MMU)
*************************************************************/

// Remove small patches (noise)
var MMU_pixels = 50;

var cc = mangroveBuffered
  .connectedPixelCount(1024, false);

var mangroveClean = mangroveBuffered
  .updateMask(cc.gte(MMU_pixels));

Map.addLayer(
  mangroveClean,
  {palette: ['800000']},
  "Mangrove Buffered + MMU",
  false
);


/************************************************************
 CONVERT RASTER TO VECTOR (GENERALIZED BOUNDARY)
*************************************************************/

var mangroveVec = mangroveClean.reduceToVectors({
  geometry: aoi,
  scale: 30,              // coarser scale for lightweight geometry
  geometryType: "polygon",
  eightConnected: false,
  bestEffort: true,
  maxPixels: 1e13
});

Map.addLayer(
  mangroveVec,
  {color: "cyan"},
  "Final Mangrove Boundary (Vector)"
);


/************************************************************
 EXPORT OPTIONS
*************************************************************/

// Export to Google Drive (Shapefile)
Export.table.toDrive({
  collection: mangroveVec,
  description: "---",
  fileNamePrefix: "---",
  fileFormat: "SHP"
});

// Optional: Export to Earth Engine Asset
// Export.table.toAsset({
//   collection: mangroveVec,
//   description: "---",
//   assetId: "your asset"
// });
