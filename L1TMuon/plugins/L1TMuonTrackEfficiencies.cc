
// 
// Class: L1TMuonTrackEfficiencies
//
// Info: Processes a track into histograms of delta-phis and such
//
// Author: L. Gray (FNAL)
//

#include <memory>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitive.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitiveFwd.h"

#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrack.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"

// detector ID types
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"

using namespace L1TMuon;

typedef edm::ParameterSet PSet;

class L1TMuonTrackEfficiencies : public edm::EDAnalyzer {  
  typedef TH1F* hist1dfp;
  typedef TH2F* hist2dfp;
  typedef std::map<std::string,hist1dfp> hist1dmap;
  typedef std::map<std::string,hist2dfp> hist2dmap;
  typedef std::map<std::string,double> deltas_map;
  typedef std::vector<double> vdouble;
  
public:
  L1TMuonTrackEfficiencies(const PSet&);
  ~L1TMuonTrackEfficiencies() {}

  void analyze(const edm::Event&, const edm::EventSetup&);  
private:  
  edm::Service<TFileService> _fs;
  edm::InputTag _geninput;
  edm::InputTag _trkInput;
  vdouble thresholds;
  vdouble sim_pt_bins;
  std::vector<unsigned> allowed_modes;
  // gen histograms
  hist1dfp gen_eta, gen_phi, gen_pt;  
  // track histograms
  hist1dfp denominator;
  std::vector<hist1dfp> numerators;
  std::vector<hist1dfp> oopt_residuals;

  edm::ESHandle<L1MuTriggerScales> position_scales;
  edm::ESHandle<L1MuTriggerPtScale> pt_scale;
  unsigned long long position_scalesCacheID;
  unsigned long long pt_scaleCacheID;
};

L1TMuonTrackEfficiencies::L1TMuonTrackEfficiencies(const PSet& p) {
  _geninput = p.getParameter<edm::InputTag>("genSrc");
  _trkInput = p.getParameter<edm::InputTag>("trackSrc");

  gen_eta = _fs->make<TH1F>("hgen_eta","Generated Eta",100,-2.5,2.5);
  gen_phi = _fs->make<TH1F>("hgen_phi","Generated Eta",100,-M_PI,M_PI);
  gen_pt = _fs->make<TH1F>("hgen_pt","Generated Eta",400,0,200);
  
  vdouble binning = p.getParameter<vdouble>("binning");

  thresholds = p.getParameter<vdouble>("thresholds");
  for( auto& threshold : thresholds ) {
    hist1dfp hist = _fs->make<TH1F>(Form("numerator_thresh%.3f",threshold),
				    ";Gen p_{T} (GeV)",
				    binning.size()-1,binning.data());
    numerators.push_back(hist);    
  }
  denominator = _fs->make<TH1F>("denominator",";Gen p_{T} (GeV)",
				binning.size()-1,binning.data());  

  sim_pt_bins = p.getParameter<vdouble>("sim_pt_bins");
  sim_pt_bins.push_back(1000);
  for(unsigned i = 0; i < sim_pt_bins.size()-1; ++i ) {
    hist1dfp hist2 = _fs->make<TH1F>(Form("ooptresidual_%.3f_%.3f",
					  sim_pt_bins[i],sim_pt_bins[i+1]),
				     ";(1/p_{T}^{L1} - 1/p_{T}^{gen})/(1/p_{T}^{gen})" ,
				     10000,-100,100);
    oopt_residuals.push_back(hist2);
  }
  allowed_modes = p.getParameter<std::vector<unsigned> >("allowed_modes");
}

void L1TMuonTrackEfficiencies::analyze(const edm::Event& ev, 
				       const edm::EventSetup& es) {
  const unsigned long long poscache = 
    es.get<L1MuTriggerScalesRcd>().cacheIdentifier();
  const unsigned long long ptcache = 
    es.get<L1MuTriggerPtScaleRcd>().cacheIdentifier();
  if(  poscache != position_scalesCacheID || ptcache  != pt_scaleCacheID ) {
    es.get<L1MuTriggerScalesRcd>().get(position_scales);
    es.get<L1MuTriggerPtScaleRcd>().get(pt_scale);
    
    position_scalesCacheID = poscache;
    pt_scaleCacheID = ptcache;
  }


  //this just handles single muons for now
  edm::Handle<reco::GenParticleCollection> genps;
  ev.getByLabel(_geninput,genps);
  auto bgen = genps->cbegin();
  auto egen = genps->cend();
  for( ; bgen != egen; ++bgen ) {
    if( std::abs(bgen->pdgId()) == 13 ) {
      gen_eta->Fill(bgen->eta());
      gen_phi->Fill(bgen->phi());
      gen_pt->Fill(bgen->pt());
      
      edm::Handle<InternalTrackCollection> trks;
      ev.getByLabel(_trkInput,trks);
      
      InternalTrackRef best_track;
      int best_pt = -1;
      auto tbeg = trks->cbegin();
      auto trk  = tbeg;
      auto tend = trks->cend();

      for( ; trk != tend; ++trk ){	
	if( (int)trk->pt_packed() > best_pt ) {
	  best_track = InternalTrackRef(trks,trk-tbeg);
	  best_pt = best_track->pt_packed();
	}
      }

      if( best_track.isNonnull() && std::abs(best_track->wheel()) <= 4) {
	denominator->Fill(bgen->pt());
	std::vector<unsigned>::const_iterator has_mode = 
	  std::find(allowed_modes.begin(),
		    allowed_modes.end(),
		    best_track->dtMode());
	std::vector<unsigned>::const_iterator has_cscmode = std::find(allowed_modes.begin(),
								      allowed_modes.end(),
								      best_track->cscMode());
	
	if( (allowed_modes.size() || has_mode != allowed_modes.end() || has_cscmode != allowed_modes.end()) &&
	    best_track->quality() <= 3 && best_track->quality() > 0    ) {
	  int pt_packed = best_track->pt_packed();
	  double ptValue = pt_scale->getPtScale()->getCenter(pt_packed);
	  for( unsigned i = 0; i < thresholds.size(); ++i ) {
	    if( thresholds[i] < ptValue ) {
	      numerators[i]->Fill(bgen->pt());	    
	    }
	  }	
	  for( unsigned i = 0; i < sim_pt_bins.size() - 1; ++i ) {
	    if( bgen->pt() > sim_pt_bins[i] && 
		bgen->pt() < sim_pt_bins[i+1] ) {
	      oopt_residuals[i]->Fill(bgen->pt()/ptValue - 1.0);
	    }
	  }
	}
      }    
    }
  }  
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonTrackEfficiencies);
