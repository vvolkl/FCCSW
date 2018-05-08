#ifndef RECTRACKER_FCCHHSEEDINGGRAPHTOOL_H
#define RECTRACKER_FCCHHSEEDINGGRAPHTOOL_H

// from Gaudi
#include "GaudiAlg/GaudiTool.h"

// FCCSW
#include "FWCore/DataHandle.h"
#include "RecInterface/ILayerGraphTool.h"

#include "tricktrack/CMGraph.h"


class FCChhSeedingGraphTool : public GaudiTool, virtual public ILayerGraphTool {
public:
  FCChhSeedingGraphTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~FCChhSeedingGraphTool() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode finalize() override final;
  virtual tricktrack::CMGraph graph() override final;
};

#endif /* RECTRACKER_FCCHHSEEDINGGRAPHTOOL_H */
