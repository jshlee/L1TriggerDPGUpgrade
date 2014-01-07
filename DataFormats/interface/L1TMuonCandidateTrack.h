#ifndef __L1TMUON_CANDIDATETRACK_H__
#define __L1TMUON_CANDIDATETRACK_H__
// 
// Class: L1TMuon::CandidateTrack
//
// Info: This class represents (one of the) final tracks output by
//       L1ITMu after sorting. It is just a L1MuGMTCand with a few
//       extra frills.
//
// Author: L. Gray (FNAL)
//

#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitiveFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitive.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrackFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrack.h"
#include "DataFormats/Common/interface/Ref.h"

namespace L1TMuon{
  
  class CandidateTrack : public L1MuGMTCand {
  public:
    CandidateTrack() {}
    ~CandidateTrack() {}

    CandidateTrack(const InternalTrackRef&);
    
    InternalTrackRef parent() const { return _parent; }
             
    const TriggerPrimitiveStationMap& getStubs() const 
      { return _parent->getStubs(); }

    unsigned long mode()     const { return _parent->mode(); }
    unsigned long dtMode()   const { return _parent->dtMode(); }
    unsigned long cscMode()  const { return _parent->cscMode(); }
    unsigned long rpcbMode() const { return _parent->rpcbMode(); }
    unsigned long rpcfMode() const { return _parent->rpcfMode(); }
    unsigned long gemMode()  const { return _parent->gemMode(); }

  private:    
    InternalTrackRef _parent;
  };
}

#endif
