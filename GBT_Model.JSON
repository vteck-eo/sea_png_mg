/**************************************************************
 🌿 Mangrove Mapping using Satellite Embeddings in GEE
 Description:
 Embedding-based mangrove classification using
 GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL
**************************************************************/


/**************************************************************
 DEFINE COASTAL AOI
 - Shoreline buffer
 - Intersection with national boundary
**************************************************************/

// Load shoreline dataset
var shoreline = ee.FeatureCollection("your_AOI_asset");

// Filter shoreline using attribute
var shorelineFiltered = shoreline.filter(
  ee.Filter.eq("Your_Attribute", "Your_Object")
);

Map.addLayer(shorelineFiltered, {color:"blue"}, "Shoreline", false);

// Create coastal buffer (30 km)
var bufferDist = 30000;  // meters

var shorelineBuffer = shorelineFiltered
  .map(function(f){ return f.buffer(bufferDist); })
  .union();

Map.addLayer(
  shorelineBuffer,
  {color:"red"},
  "Shoreline Buffer (30 km)",
  false
);

// Load FAO boundary (GAUL 2015)
var fao = ee.FeatureCollection("FAO/GAUL/2015/level0");

var countryBoundary = fao
  .filter(ee.Filter.eq('ADM0_NAME', 'Country'))
  .geometry();

// Final AOI = Country ∩ Coastal Buffer
var aoi = countryBoundary.intersection(shorelineBuffer, 1);

Map.addLayer(aoi, {color: 'cyan'}, 'Final AOI', false);


/**************************************************************
 LOAD TRAINING SAMPLES
**************************************************************/

var samples = ee.FeatureCollection("your_sample_asset")
  .filterBounds(aoi)
  .filter(ee.Filter.neq("A01", null));  // remove invalid samples

print('Sample preview:', samples.limit(10));
print('Class histogram:', samples.aggregate_histogram('land_class'));


/**************************************************************
 LOAD GOOGLE SATELLITE EMBEDDINGS (ANNUAL)
**************************************************************/

var annualEmbeddings =
  ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL');

function getYearlyEmbeddings(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');

  return annualEmbeddings
    .filterDate(startDate, endDate)
    .filterBounds(aoi)
    .mosaic()
    .clip(aoi);
}

// Example year
var img = getYearlyEmbeddings(2017);
print("Embedding image:", img);


/**************************************************************
 PREPARE BALANCED TRAINING DATA
**************************************************************/

var mangrove = samples
  .filter(ee.Filter.eq('land_class', 1))
  .randomColumn('rand')
  .sort('rand')
  .limit(1000);

var nonMangrove = samples
  .filter(ee.Filter.eq('land_class', 0))
  .randomColumn('rand')
  .sort('rand')
  .limit(1000);

var balanced = mangrove.merge(nonMangrove);

Map.addLayer(balanced, {color: 'yellow'}, 'Balanced Samples', false);


/**************************************************************
 EXTRACT EMBEDDING FEATURES
**************************************************************/

var trainingData = img.sampleRegions({
  collection: balanced,
  properties: ['land_class'],
  scale: 10,
  geometries: true
}).filter(ee.Filter.notNull(img.bandNames()));


/**************************************************************
 TRAIN GRADIENT TREE BOOST CLASSIFIER
**************************************************************/

var withRandom = trainingData.randomColumn('rand');
var split = 0.95;

var training = withRandom.filter(ee.Filter.lt('rand', split));
var testing  = withRandom.filter(ee.Filter.gte('rand', split));

var model = ee.Classifier.smileGradientTreeBoost({
  numberOfTrees: 100,
  shrinkage: 0.075
}).train({
  features: training,
  classProperty: 'land_class',
  inputProperties: img.bandNames()
}).setOutputMode('PROBABILITY');

print("Model explanation:", model.explain());


/**************************************************************
 FEATURE IMPORTANCE
**************************************************************/

var importance = ee.Dictionary(model.explain().get('importance'));

var importanceFC = ee.FeatureCollection(
  importance.keys().map(function(k){
    return ee.Feature(null, {
      band: k,
      importance: importance.get(k)
    });
  })
);

var sorted = importanceFC.sort('importance', false);

var chart = ui.Chart.feature.byFeature(
  sorted, 'band', 'importance'
).setChartType('ColumnChart')
 .setOptions({
   title: 'Embedding Band Importance',
   legend: {position: 'none'}
 });

print(chart);


/**************************************************************
 GENERATE MANGROVE PROBABILITY MAP
**************************************************************/

var probability = img.classify(model).multiply(100);

Map.addLayer(probability, {
  min: 0,
  max: 100,
  palette: ['white', 'yellow', 'orange', 'red', 'purple']
}, 'Mangrove Probability');


/**************************************************************
 ACCURACY ASSESSMENT
**************************************************************/

var tested = testing.classify(model);

var validated = tested.map(function(f) {
  var prob = ee.Number(f.get('classification'));
  var pred = prob.gt(0.5).int();  // threshold = 0.5
  return f.set('predicted', pred);
});

var cm = validated.errorMatrix('land_class', 'predicted');

print('Confusion Matrix:', cm);
print('Overall Accuracy:', cm.accuracy());
print('Producer Accuracy:', cm.producersAccuracy());
print('User Accuracy:', cm.consumersAccuracy());
