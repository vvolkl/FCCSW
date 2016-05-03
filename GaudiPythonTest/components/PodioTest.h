
#include "GaudiAlg/GaudiAlgorithm.h"

#include "FWCore/DataHandle.h"

// datamodel
namespace fcc {
class TrackClusterCollection;
class TrackCluster;
}

class PodioTest: public GaudiAlgorithm {
public:
  explicit PodioTest(const std::string& ,ISvcLocator*);
  virtual ~PodioTest();
  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();
private:
  DataHandle<fcc::TrackClusterCollection> m_trackClusters;
};

