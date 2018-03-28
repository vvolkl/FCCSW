#ifndef RECINTERFACE_ITRACKFITTINGTOOL_H
#define RECINTERFACE_ITRACKFITTINGTOOL_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"

#include <map>

namespace fcc {
class PositionedTrackHitCollection;
class TrackCollection;
class TrackStateCollection;
}


class ITrackFittingTool : virtual public IAlgTool {
public:
  DeclareInterfaceID(ITrackFittingTool, 1, 0);
  virtual std::pair<fcc::TrackCollection*, fcc::TrackStateCollection*> fitTracks(const fcc::PositionedTrackHitCollection* theHits, std::multimap<unsigned int, unsigned int> seedmap) = 0;

};

#endif /* RECINTERFACE_ITRACKFITTINGTOOL_H */
