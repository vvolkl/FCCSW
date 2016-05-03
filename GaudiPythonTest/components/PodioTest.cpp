
// datamodel
#include "datamodel/TrackClusterCollection.h"

#include "PodioTest.h"

DECLARE_ALGORITHM_FACTORY(PodioTest)

 PodioTest::PodioTest(const std::string& aName ,ISvcLocator* aSvcLocator):
  GaudiAlgorithm(aName, aSvcLocator) { 
  declareOutput("trackClusters", m_trackClusters,"hits/trackerClusters");
}

PodioTest::~PodioTest() {}

StatusCode PodioTest::initialize() {
  return GaudiAlgorithm::initialize();
}

StatusCode PodioTest::finalize() {
  return GaudiAlgorithm::finalize();
}

StatusCode PodioTest::execute() {
  fcc::TrackClusterCollection* edmClusters = new fcc::TrackClusterCollection();
  for (int iter_coll=0; iter_coll < 100000; ++iter_coll) {
    fcc::TrackCluster edmCluster = edmClusters->create();
    fcc::BareCluster& edmClusterCore = edmCluster.Core();
    // fill with some values for testing
    edmClusterCore.position.X = iter_coll; 
    edmClusterCore.position.Y = 2; 
    edmClusterCore.position.Z = 3; 
    edmClusterCore.Energy = 4;
    edmClusterCore.Time = 5;
   }
  m_trackClusters.put(edmClusters);
  return StatusCode::SUCCESS;
}
