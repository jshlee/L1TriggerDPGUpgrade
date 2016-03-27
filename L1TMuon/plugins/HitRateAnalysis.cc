// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMCSCPadDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMCSCCoPadDigiCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TString.h"

using namespace std;
using namespace edm;

class HitRateAnalysis : public edm::EDAnalyzer {
public:
  explicit HitRateAnalysis(const edm::ParameterSet&);
  ~HitRateAnalysis();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag gemDigiInput_;
  edm::InputTag gemPadDigiInput_;
  edm::InputTag gemCoPadDigiInput_;

  edm::Service<TFileService> fs;

  TH1F* h_digiPerNBX;
  TH1F* h_padPerNBX;
  TH1F* h_coPadPerNBX;
  int m_bxnumber;
  const int m_maxbx = 40;
  
  int m_digiPerNBX[73][2][10][192];
  int m_padPerNBX[73][2][10][192];
  int m_coPadPerNBX[73][2][10][192];

  TH1F* h_nchamber;
  
};
HitRateAnalysis::HitRateAnalysis(const edm::ParameterSet& iConfig)
{
  gemDigiInput_ = iConfig.getParameter<edm::InputTag>("gemDigiInput");
  gemPadDigiInput_ = iConfig.getParameter<edm::InputTag>("gemPadDigiInput");
  gemCoPadDigiInput_ = iConfig.getParameter<edm::InputTag>("gemCoPadDigiInput");

  m_bxnumber = 0;
  
  h_nchamber=fs->make<TH1F>("nchamber","nchamber",73,0,73);
  h_nchamber->GetXaxis()->SetTitle(Form("same digi hits within %2i BX", m_maxbx));
  h_nchamber->GetYaxis()->SetTitle("Counts");

  h_digiPerNBX=fs->make<TH1F>("digiPerNBX","digiPerNBX",10,0,10);
  h_digiPerNBX->GetXaxis()->SetTitle(Form("same digi hits within %2i BX", m_maxbx));
  h_digiPerNBX->GetYaxis()->SetTitle("Counts");
  
  h_padPerNBX=fs->make<TH1F>("padPerNBX","padPerNBX",10,0,10);
  h_padPerNBX->GetXaxis()->SetTitle(Form("same pad hits within %2i BX", m_maxbx));
  h_padPerNBX->GetYaxis()->SetTitle("Counts");
  
  h_coPadPerNBX=fs->make<TH1F>("coPadPerNBX","coPadPerNBX",10,0,10);
  h_coPadPerNBX->GetXaxis()->SetTitle(Form("same coPad hits within %2i BX", m_maxbx));
  h_coPadPerNBX->GetYaxis()->SetTitle("Counts");

}
HitRateAnalysis::~HitRateAnalysis(){}
void
HitRateAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GEMDigiCollection> gem_digis;  
  edm::Handle<GEMCSCPadDigiCollection> gemcscpad_digis;
  edm::Handle<GEMCSCCoPadDigiCollection> gemcsccopad_digis;
  iEvent.getByLabel(gemDigiInput_, gem_digis);
  iEvent.getByLabel(gemPadDigiInput_, gemcscpad_digis);
  iEvent.getByLabel(gemCoPadDigiInput_, gemcsccopad_digis);

  if (m_bxnumber == 0){
    for (int i = 0; i < 73; ++i)
      for (int j = 0; j < 2; ++j)
	for (int f = 0; f < 10; ++f)
	  for (int k = 0; k < 192; ++k){
	    m_digiPerNBX[i][j][f][k] = 0;
	    m_padPerNBX[i][j][f][k] = 0;
	    m_coPadPerNBX[i][j][f][k] = 0;
	  }
  }
  
  for(GEMDigiCollection::DigiRangeIterator cItr = gem_digis->begin(); cItr != gem_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1){
	m_digiPerNBX[36+id.region()*id.chamber()][id.layer()-1][id.roll()-1][digiIt->strip()]++;
	h_nchamber->Fill(id.roll()-1);
      }
    }
  }
  
  for(GEMCSCPadDigiCollection::DigiRangeIterator cItr = gemcscpad_digis->begin(); cItr != gemcscpad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMCSCPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1){	
	m_padPerNBX[36+id.region()*id.chamber()][id.layer()-1][id.roll()-1][digiIt->pad()]++;
      }
    }
  }

  for(GEMCSCCoPadDigiCollection::DigiRangeIterator cItr = gemcsccopad_digis->begin(); cItr != gemcsccopad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMCSCCoPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1){
	m_coPadPerNBX[36+id.region()*id.chamber()][id.layer()-1][id.roll()-1][digiIt->first().pad()]++;
	//cout <<"id.layer() " << id.layer()-1 << endl;
      }
    }
  }
  m_bxnumber++;
  
  if (m_bxnumber == m_maxbx){
    for (int i = 0; i < 73; ++i)
      if (i != 36)
	for (int j = 0; j < 2; ++j)
	  for (int f = 0; f < 10; ++f)
	    for (int k = 0; k < 192; ++k){
	      h_digiPerNBX->Fill(m_digiPerNBX[i][j][f][k]);
	      h_padPerNBX->Fill(m_padPerNBX[i][j][f][k]);
	      if (j==0)
		h_coPadPerNBX->Fill(m_coPadPerNBX[i][j][f][k]);
	    }
    m_bxnumber = 0;
  }
  
}

void HitRateAnalysis::beginJob(){}
void HitRateAnalysis::endJob(){}

//define this as a plug-in
DEFINE_FWK_MODULE(HitRateAnalysis);
