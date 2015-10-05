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

using namespace std;
using namespace edm;

class HitAnalysis : public edm::EDAnalyzer {
public:
  explicit HitAnalysis(const edm::ParameterSet&);
  ~HitAnalysis();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag gemDigiInput_;
  edm::InputTag gemPadDigiInput_;
  edm::InputTag gemCoPadDigiInput_;

  edm::Service<TFileService> fs;

  TH1F* h_nPadPerChamber;
  TH2F* h_nPadPerChamberPerBX;
  TH1F* h_nPadPerChamber10;
  TH2F* h_nPadPerChamberPerBX10;
  TH1F* h_nPadPerChamberGE11;
  TH2F* h_nPadPerChamberPerBXGE11;

  TH1F* h_nPadChamber;
  TH2F* h_nPadChamberPerBX;
  TH1F* h_nPadChamber10;
  TH2F* h_nPadChamberPerBX10;
  TH1F* h_nPadChamberGE11;
  TH2F* h_nPadChamberPerBXGE11;

  TH1F* h_nPadRangeChamber;
  TH2F* h_nPadRangeChamberPerBX;
  TH1F* h_nPadRangeChamber10;
  TH2F* h_nPadRangeChamberPerBX10;
  TH1F* h_nPadRangeChamberGE11;
  TH2F* h_nPadRangeChamberPerBXGE11;
  
  TH1F* h_nCoPadPerChamber;
  TH2F* h_nCoPadPerChamberPerBX;
  TH1F* h_nCoPadPerChamber10;
  TH2F* h_nCoPadPerChamberPerBX10;
  TH1F* h_nCoPadPerChamberGE11;
  TH2F* h_nCoPadPerChamberPerBXGE11;

  TH2F* h_maxPadvsCoPad;
  TH2F* h_maxPadvsCoPad10;
  TH2F* h_maxPadvsCoPadGE11;
  
};
HitAnalysis::HitAnalysis(const edm::ParameterSet& iConfig)
{
  gemDigiInput_ = iConfig.getParameter<edm::InputTag>("gemDigiInput");
  gemPadDigiInput_ = iConfig.getParameter<edm::InputTag>("gemPadDigiInput");
  gemCoPadDigiInput_ = iConfig.getParameter<edm::InputTag>("gemCoPadDigiInput");

  h_nPadPerChamber=fs->make<TH1F>("nPadPerChamber20","nPadPerChamber20",100,0,100);
  h_nPadPerChamber->GetXaxis()->SetTitle("no. Pads Per Chamber");
  h_nPadPerChamber->GetYaxis()->SetTitle("Counts");
  h_nPadPerChamberPerBX=fs->make<TH2F>("nPadPerChamberPerBX20","nPadPerChamberPerBX20",50,-25,25,100,0,100);
  h_nPadPerChamberPerBX->GetXaxis()->SetTitle("BX");
  h_nPadPerChamberPerBX->GetYaxis()->SetTitle("no. Pads Per Chamber");
  h_nPadPerChamber10=fs->make<TH1F>("nPadPerChamber10","nPadPerChamber10",100,0,100);
  h_nPadPerChamber10->GetXaxis()->SetTitle("no. Pads Per Chamber");
  h_nPadPerChamber10->GetYaxis()->SetTitle("Counts");
  h_nPadPerChamberPerBX10=fs->make<TH2F>("nPadPerChamberPerBX10","nPadPerChamberPerBX10",50,-25,25,100,0,100);
  h_nPadPerChamberPerBX10->GetXaxis()->SetTitle("BX");
  h_nPadPerChamberPerBX10->GetYaxis()->SetTitle("no. Pads Per Chamber");
  h_nPadPerChamberGE11=fs->make<TH1F>("nPadPerChamberGE11","nPadPerChamberGE11",100,0,100);
  h_nPadPerChamberGE11->GetXaxis()->SetTitle("no. Pads Per Chamber");
  h_nPadPerChamberGE11->GetYaxis()->SetTitle("Counts");
  h_nPadPerChamberPerBXGE11=fs->make<TH2F>("nPadPerChamberPerBXGE11","nPadPerChamberPerBXGE11",50,-25,25,100,0,100);
  h_nPadPerChamberPerBXGE11->GetXaxis()->SetTitle("BX");
  h_nPadPerChamberPerBXGE11->GetYaxis()->SetTitle("no. Pads Per Chamber");

  h_nCoPadPerChamber=fs->make<TH1F>("nCoPadPerChamber20","nCoPadPerChamber20",100,0,100);
  h_nCoPadPerChamber->GetXaxis()->SetTitle("no. CoPads Per Chamber");
  h_nCoPadPerChamber->GetYaxis()->SetTitle("Counts");
  h_nCoPadPerChamberPerBX=fs->make<TH2F>("nCoPadPerChamberPerBX20","nCoPadPerChamberPerBX20",50,-25,25,100,0,100);
  h_nCoPadPerChamberPerBX->GetXaxis()->SetTitle("BX");
  h_nCoPadPerChamberPerBX->GetYaxis()->SetTitle("no. CoPads Per Chamber");
  h_nCoPadPerChamber10=fs->make<TH1F>("nCoPadPerChamber10","nCoPadPerChamber10",100,0,100);
  h_nCoPadPerChamber10->GetXaxis()->SetTitle("no. CoPads Per Chamber");
  h_nCoPadPerChamber10->GetYaxis()->SetTitle("Counts");
  h_nCoPadPerChamberPerBX10=fs->make<TH2F>("nCoPadPerChamberPerBX10","nCoPadPerChamberPerBX10",50,-25,25,100,0,100);
  h_nCoPadPerChamberPerBX10->GetXaxis()->SetTitle("BX");
  h_nCoPadPerChamberPerBX10->GetYaxis()->SetTitle("no. CoPads Per Chamber");
  h_nCoPadPerChamberGE11=fs->make<TH1F>("nCoPadPerChamberGE11","nCoPadPerChamberGE11",100,0,100);
  h_nCoPadPerChamberGE11->GetXaxis()->SetTitle("no. CoPads Per Chamber");
  h_nCoPadPerChamberGE11->GetYaxis()->SetTitle("Counts");
  h_nCoPadPerChamberPerBXGE11=fs->make<TH2F>("nCoPadPerChamberPerBXGE11","nCoPadPerChamberPerBXGE11",50,-25,25,100,0,100);
  h_nCoPadPerChamberPerBXGE11->GetXaxis()->SetTitle("BX");
  h_nCoPadPerChamberPerBXGE11->GetYaxis()->SetTitle("no. CoPads Per Chamber");

  h_nPadChamber=fs->make<TH1F>("nPadChamber20","nPadChamber20",100,0,100);
  h_nPadChamber->GetXaxis()->SetTitle("no. Pad Ranges per Chamber");
  h_nPadChamber->GetYaxis()->SetTitle("Counts");
  h_nPadChamberPerBX=fs->make<TH2F>("nPadChamberPerBX20","nPadChamberPerBX20",50,-25,25,100,0,100);
  h_nPadChamberPerBX->GetXaxis()->SetTitle("BX");
  h_nPadChamberPerBX->GetYaxis()->SetTitle("no. Pad Ranges per Chamber");
  h_nPadChamber10=fs->make<TH1F>("nPadChamber10","nPadChamber10",100,0,100);
  h_nPadChamber10->GetXaxis()->SetTitle("no. Pad Ranges per Chamber");
  h_nPadChamber10->GetYaxis()->SetTitle("Counts");
  h_nPadChamberPerBX10=fs->make<TH2F>("nPadChamberPerBX10","nPadChamberPerBX10",50,-25,25,100,0,100);
  h_nPadChamberPerBX10->GetXaxis()->SetTitle("BX");
  h_nPadChamberPerBX10->GetYaxis()->SetTitle("no. Pad Ranges per Chamber");
  h_nPadChamberGE11=fs->make<TH1F>("nPadChamberGE11","nPadChamberGE11",100,0,100);
  h_nPadChamberGE11->GetXaxis()->SetTitle("no. Pad Ranges per Chamber");
  h_nPadChamberGE11->GetYaxis()->SetTitle("Counts");
  h_nPadChamberPerBXGE11=fs->make<TH2F>("nPadChamberPerBXGE11","nPadChamberPerBXGE11",50,-25,25,100,0,100);
  h_nPadChamberPerBXGE11->GetXaxis()->SetTitle("BX");
  h_nPadChamberPerBXGE11->GetYaxis()->SetTitle("no. Pad Ranges per Chamber");

  h_nPadRangeChamber=fs->make<TH1F>("nPadRangeChamber20","nPadRangeChamber20",100,0,100);
  h_nPadRangeChamber->GetXaxis()->SetTitle("size of Pad Range");
  h_nPadRangeChamber->GetYaxis()->SetTitle("Counts");
  h_nPadRangeChamberPerBX=fs->make<TH2F>("nPadRangeChamberPerBX20","nPadRangeChamberPerBX20",50,-25,25,100,0,100);
  h_nPadRangeChamberPerBX->GetXaxis()->SetTitle("BX");
  h_nPadRangeChamberPerBX->GetYaxis()->SetTitle("size of Pad Range");
  h_nPadRangeChamber10=fs->make<TH1F>("nPadRangeChamber10","nPadRangeChamber10",100,0,100);
  h_nPadRangeChamber10->GetXaxis()->SetTitle("size of Pad Range");
  h_nPadRangeChamber10->GetYaxis()->SetTitle("Counts");
  h_nPadRangeChamberPerBX10=fs->make<TH2F>("nPadRangeChamberPerBX10","nPadRangeChamberPerBX10",50,-25,25,100,0,100);
  h_nPadRangeChamberPerBX10->GetXaxis()->SetTitle("BX");
  h_nPadRangeChamberPerBX10->GetYaxis()->SetTitle("size of Pad Range");
  h_nPadRangeChamberGE11=fs->make<TH1F>("nPadRangeChamberGE11","nPadRangeChamberGE11",100,0,100);
  h_nPadRangeChamberGE11->GetXaxis()->SetTitle("size of Pad Range");
  h_nPadRangeChamberGE11->GetYaxis()->SetTitle("Counts");
  h_nPadRangeChamberPerBXGE11=fs->make<TH2F>("nPadRangeChamberPerBXGE11","nPadRangeChamberPerBXGE11",50,-25,25,100,0,100);
  h_nPadRangeChamberPerBXGE11->GetXaxis()->SetTitle("BX");
  h_nPadRangeChamberPerBXGE11->GetYaxis()->SetTitle("size of Pad Range");

  h_maxPadvsCoPad=fs->make<TH2F>("maxPadvsCoPad","maxPadvsCoPad",100,0,100,100,0,100);
  h_maxPadvsCoPad->GetXaxis()->SetTitle("max no. Pads");
  h_maxPadvsCoPad->GetYaxis()->SetTitle("no. CoPads");

  h_maxPadvsCoPad10=fs->make<TH2F>("maxPadvsCoPad10","maxPadvsCoPad10",100,0,100,100,0,100);
  h_maxPadvsCoPad10->GetXaxis()->SetTitle("max no. Pads");
  h_maxPadvsCoPad10->GetYaxis()->SetTitle("no. CoPads");
  
  h_maxPadvsCoPadGE11=fs->make<TH2F>("maxPadvsCoPadGE11","maxPadvsCoPadGE11",100,0,100,100,0,100);
  h_maxPadvsCoPadGE11->GetXaxis()->SetTitle("max no. Pads");
  h_maxPadvsCoPadGE11->GetYaxis()->SetTitle("no. CoPads");

}
HitAnalysis::~HitAnalysis(){}
void
HitAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GEMDigiCollection> gem_digis;  
  edm::Handle<GEMCSCPadDigiCollection> gemcscpad_digis;
  edm::Handle<GEMCSCCoPadDigiCollection> gemcsccopad_digis;
  iEvent.getByLabel(gemDigiInput_, gem_digis);
  iEvent.getByLabel(gemPadDigiInput_, gemcscpad_digis);
  iEvent.getByLabel(gemCoPadDigiInput_, gemcsccopad_digis);

  int nPadPerChamber[73][2][40] = {};
  int nPadPerChamber10[73][2][40] = {};
  int nPadPerChamberGE11[73][2][40] = {};

  int nPadChamber[73][2][192] = {};
  int nPadChamber10deg[73][2][192] = {};
  int nPadChamberGE11deg[73][2][192] = {};

  for(GEMCSCPadDigiCollection::DigiRangeIterator cItr = gemcscpad_digis->begin(); cItr != gemcscpad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);

    for (GEMCSCPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1){
	nPadPerChamberGE11[36+id.region()*id.chamber()][id.layer()-1][15+digiIt->bx()]++;
	if (digiIt->bx() == 0) nPadChamberGE11deg[36+id.region()*id.chamber()][id.layer()-1][digiIt->pad()]++;
	//cout << "ge11 id.layer()" << id.layer() << endl;
      }
      if (id.station() == 2){
	//cout << "ge21 id.layer()" << id.layer() << endl;
	nPadPerChamber[36+id.region()*id.chamber()][id.layer()-1][15+digiIt->bx()]++;
	if (digiIt->bx() == 0) nPadChamber[36+id.region()*id.chamber()][id.layer()-1][digiIt->pad()]++;
	if (digiIt->pad() > 96){
	  nPadPerChamber10[36+id.region()*id.chamber()][id.layer()-1][15+digiIt->bx()]++;
	  if (digiIt->bx() == 0) nPadChamber10deg[36+id.region()*id.chamber()][id.layer()-1][digiIt->pad()]++;
	}
      }
    }
  }
  for (int i = 0; i < 73; i++){
    for (int h = 0; h < 2; h++){
      
      if ( i == 36 ) continue;
      //if (nPadPerChamber[i][h][15])
      if ( i >= 18 && i <= 54)
	h_nPadPerChamber->Fill(nPadPerChamber[i][h][15]);
      //    if (nPadPerChamber10[i][h][15])
      if ( i >= 18 && i <= 54)
	h_nPadPerChamber10->Fill(nPadPerChamber10[i][h][15]);
      //    if (nPadPerChamberGE11[i][h][15])
      h_nPadPerChamberGE11->Fill(nPadPerChamberGE11[i][h][15]);
      for (int j = -15; j < 25; j++){
	//      if (nPadPerChamber[i][h][15+j])
	if ( i >= 18 && i <= 54)
	  h_nPadPerChamberPerBX->Fill(j,nPadPerChamber[i][h][15+j]);
	//      if (nPadPerChamber10[i][h][15+j])
	if ( i >= 18 && i <= 54)
	  h_nPadPerChamberPerBX10->Fill(j,nPadPerChamber10[i][h][15+j]);
	//      if (nPadPerChamberGE11[i][h][15+j])
	h_nPadPerChamberPerBXGE11->Fill(j,nPadPerChamberGE11[i][h][15+j]);
      }
      int nContPadRanges = 0;
      int nContPadRanges10 = 0;
      int nContPadRangesGE11 = 0;
      int maxrange = 0;
      int maxrange10 = 0;
      int maxrangeGE11 = 0;
      int j=0;
      for (int k = 0; k < 192; k++){
	
	if (nPadChamber[i][h][k] && maxrange < 8){
	  maxrange++;
	}
	else if (maxrange){
	  if (j == 0) h_nPadRangeChamber->Fill(maxrange);
	  h_nPadRangeChamberPerBX->Fill(j,maxrange);
	  maxrange = 0;
	  nContPadRanges++;
	}

	if (nPadChamber10deg[i][h][k] && maxrange10 < 8){
	  maxrange10++;
	}
	else if (maxrange10){
	  if (j == 0) h_nPadRangeChamber10->Fill(maxrange10);
	  h_nPadRangeChamberPerBX10->Fill(j,maxrange10);
	  maxrange10 = 0;
	  nContPadRanges10++;
	}

	if (nPadChamberGE11deg[i][h][k] && maxrangeGE11 < 8){
	  maxrangeGE11++;
	}
	else if (maxrangeGE11){
	  if (j == 0) h_nPadRangeChamberGE11->Fill(maxrangeGE11);
	  h_nPadRangeChamberPerBXGE11->Fill(j,maxrangeGE11);
	  maxrangeGE11 = 0;
	  nContPadRangesGE11++;
	}
	
      }
      
      if ( i >= 18 && i <= 54){
	if (j == 0) h_nPadChamber->Fill(nContPadRanges);
	//h_nPadChamberPerBX->Fill(j,nContPadRanges);

	if (j == 0) h_nPadChamber10->Fill(nContPadRanges10);
	//h_nPadChamberPerBX10->Fill(j,nContPadRanges10);
      }
      if (j == 0) h_nPadChamberGE11->Fill(nContPadRangesGE11);
      //h_nPadChamberPerBXGE11->Fill(j,nContPadRangesGE11);
    }
  }
  
  int nCoPadPerChamber[73][40] = {};
  int nCoPadPerChamber10[73][40] = {};
  int nCoPadPerChamberGE11[73][40] = {};
  for(GEMCSCCoPadDigiCollection::DigiRangeIterator cItr = gemcsccopad_digis->begin(); cItr != gemcsccopad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMCSCCoPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1)
	nCoPadPerChamberGE11[36+id.region()*id.chamber()][15+digiIt->first().bx()]++;
      if (id.station() == 2){
	nCoPadPerChamber[36+id.region()*id.chamber()][15+digiIt->first().bx()]++;
	if (digiIt->first().pad() > 96 && digiIt->second().pad() > 96)
	  nCoPadPerChamber10[36+id.region()*id.chamber()][15+digiIt->first().bx()]++;
      }
    }
  }
  for (int i = 0; i < 73; i++){
    if ( i == 36 ) continue;    
    //    if (nCoPadPerChamber[i][15])
    if ( i >= 18 && i <= 54)
      h_nCoPadPerChamber->Fill(nCoPadPerChamber[i][15]);
    //    if (nCoPadPerChamber10[i][15])
    if ( i >= 18 && i <= 54)
      h_nCoPadPerChamber10->Fill(nCoPadPerChamber10[i][15]);
    //    if (nCoPadPerChamberGE11[i][15])
    h_nCoPadPerChamberGE11->Fill(nCoPadPerChamberGE11[i][15]);
    for (int j = -15; j < 25; j++){
      //      if (nCoPadPerChamber[i][15+j])
      if ( i >= 18 && i <= 54)
	h_nCoPadPerChamberPerBX->Fill(j,nCoPadPerChamber[i][15+j]);
      //      if (nCoPadPerChamber10[i][15+j])
      if ( i >= 18 && i <= 54)
	h_nCoPadPerChamberPerBX10->Fill(j,nCoPadPerChamber10[i][15+j]);
      //      if (nCoPadPerChamberGE11[i][15+j])
      h_nCoPadPerChamberPerBXGE11->Fill(j,nCoPadPerChamberGE11[i][15+j]);
    }

    if (nPadPerChamber[i][0][15] || nPadPerChamber[i][1][15])
      h_maxPadvsCoPad->Fill(std::max(nPadPerChamber[i][0][15], nPadPerChamber[i][1][15]) , nCoPadPerChamber[i][15] );

    if (nPadPerChamber10[i][0][15] || nPadPerChamber10[i][1][15])
      h_maxPadvsCoPad10->Fill(std::max(nPadPerChamber10[i][0][15], nPadPerChamber10[i][1][15]) , nCoPadPerChamber10[i][15] );

    if (nPadPerChamberGE11[i][0][15] || nPadPerChamberGE11[i][1][15])
      h_maxPadvsCoPadGE11->Fill(std::max(nPadPerChamberGE11[i][0][15], nPadPerChamberGE11[i][1][15]) , nCoPadPerChamberGE11[i][15] );
    
  }
  
}

void HitAnalysis::beginJob(){}
void HitAnalysis::endJob(){}

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalysis);
