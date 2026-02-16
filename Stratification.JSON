/************************************************************
 STEP 0 — Load Probability Images (Google Satellite Embedding)
*************************************************************/
var assetList = {
  2017: "image/2017",
  2024: "image/2024"

};


function loadProb(year) {
  return ee.Image(assetList[year])
    .select("classification")
    .divide(100)
    .rename("p" + year);
}

var p_t1 = loadProb(2017);
var p_t2 = loadProb(2024);

// QA display
Map.addLayer(
  p_t2,
  {min: 0, max: 1, palette: ['blue', 'yellow', 'red']},
  'Probability 2024',
  false
);

/************************************************************
 STEP 1 — Coastal AOI (MERGE ALL FEATURES)
*************************************************************/
var coastalFC = ee.FeatureCollection(
  "your boundary" //Derived from mangrove forest boundary  
  
);

// 🔑 Merge ALL features, then simplify ONCE (100 m)
var coastalGeom = coastalFC.geometry().simplify(100);
// QA display
//Map.centerObject(coastalGeom, 7);
Map.addLayer(coastalGeom, {color: 'white'}, 'Merged Coastal AOI', true);

// Raster mask
var coastalMask = ee.Image(1).clip(coastalGeom).selfMask();

/************************************************************
 STEP 2 — Define Strata (RAW, 10 m)
*************************************************************/
var L = 0.45;
var H = 0.55;

// Stable mangrove
var S1_raw = p_t1.multiply(p_t2)
  .gte(H)
  .updateMask(coastalMask);

// Stable non-mangrove
var S2_raw = ee.Image(1).subtract(p_t1)
  .multiply(ee.Image(1).subtract(p_t2))
  .gte(H)
  .updateMask(coastalMask);

// Loss
var S3_raw = p_t1.multiply(ee.Image(1).subtract(p_t2))
  .gte(L)
  .updateMask(coastalMask);

// Gain
var S4_raw = ee.Image(1).subtract(p_t1)
  .multiply(p_t2)
  .gte(L)
  .updateMask(coastalMask);

// Uncertain
var S5_raw = S1_raw.not()
  .and(S2_raw.not())
  .and(S3_raw.not())
  .and(S4_raw.not())
  .updateMask(coastalMask);

/************************************************************
 STEP 3 — Combine Strata
*************************************************************/
var strataRaw = ee.Image(0)
  .updateMask(coastalMask)
  .where(S1_raw, 1)
  .where(S2_raw, 2)
  .where(S3_raw, 3)
  .where(S4_raw, 4)
  .where(S5_raw, 5)
  .rename("strata");

/************************************************************
 STEP 4 — APPLY MMU (50 pixels = 0.5 ha)
*************************************************************/
var MMU_pixels = 50;

// Connected components
var cc = strataRaw.connectedPixelCount(1024, false); // 4-connected

// Apply MMU once
var strataMMU = strataRaw.updateMask(cc.gte(MMU_pixels));

// QA display
Map.addLayer(
  strataMMU,
  {min: 1, max: 5, palette: ['006400','A9A9A9','0000FF','FFA500','FF00FF']},
  'Strata (MMU 0.5 ha)',
  true
);

/************************************************************
 STEP 5 — Split Clean Strata (Binary Masks)
*************************************************************/
var S1 = strataMMU.eq(1).selfMask();
var S2 = strataMMU.eq(2).selfMask();
var S3 = strataMMU.eq(3).selfMask();
var S4 = strataMMU.eq(4).selfMask();
var S5 = strataMMU.eq(5).selfMask();

/************************************************************
 STEP 6 — AREA CALCULATION (CORRECT METHOD)
 area = ee.Image.pixelArea() × binary mask
*************************************************************/
var areaImg = ee.Image.pixelArea().divide(10000); // m² → ha

function areaOf(mask, label) {
  return ee.Feature(null, {
    stratum: label,
    area_ha: areaImg.updateMask(mask).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: coastalGeom,
      scale: 10,              // native resolution
      maxPixels: 1e13,
      tileScale: 8
    }).get("area")
  });
}

var areaTable = ee.FeatureCollection([
  areaOf(S1, "S1_Stable_Mangrove"),
  areaOf(S2, "S2_Stable_NonMangrove"),
  areaOf(S3, "S3_Loss"),
  areaOf(S4, "S4_Gain"),
  areaOf(S5, "S5_Uncertain")
]);

// Safe print
//print("Area by stratum (ha)", areaTable);

/************************************************************
 STEP 7 — EXPORT AREA TABLE (CSV)
*************************************************************/
// Export.table.toDrive({
//   collection: areaTable,
//   description: "----",
//   folder: "----",
//   fileNamePrefix: "----",
//   fileFormat: "CSV",
//   selectors: ["stratum", "area_ha"]
// });



/************************************************************
  5. STRATIFIED SAMPLING (FINAL, MMU-BASED) 
*************************************************************/

// Final sample allocation (from your table)
var classValues  = [1, 2, 3, 4, 5];
var classPoints  = [119, 120, 43, 43, 43]; // Derrived from Sample allocation //https://doi.org/10.1016/j.rse.2014.02.015

var samples = strataMMU.stratifiedSample({
  numPoints: 0,                 // use classPoints
  classBand: "strata",
  region: coastalGeom,          // merged + simplified geometry
  scale: 10,                    // native resolution
  seed: 2025,
  classValues: classValues,
  classPoints: classPoints,
  geometries: true,
  dropNulls: true
});

/************************************************************
  6. ADD IDS + EXPORT
*************************************************************/

var nameDict = ee.Dictionary({
  1: "S1_Stable_Mangrove",
  2: "S2_Stable_NonMangrove",
  3: "S3_Potential_Loss",
  4: "S4_Potential_Gain",
  5: "S5_Uncertain"
});

samples = samples.map(function(f){
  var idx = ee.Number.parse(f.id()).format('%06d');
  var xy = f.geometry().coordinates();
  return f.set({
    id: idx,
    point_id: ee.String("1_").cat(idx),
    strata_name: nameDict.get(f.get("strata")),
    longitude: xy.get(0),
    latitude: xy.get(1)
  });
});

/************************************************************
  7. EXPORT SAMPLES (CSV)
*************************************************************/
Export.table.toDrive({
  collection: samples,
  description: "---",
  folder: "----",
  fileNamePrefix: "----",
  fileFormat: "CSV",
  selectors: [
    "id","point_id",
    "strata","strata_name",
    "longitude","latitude"
  ]
});
