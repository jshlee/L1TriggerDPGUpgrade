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

  typedef pair<GEMDetId,int> gemstrip;
  vector<gemstrip> m_gemstrips;
  
  edm::Service<TFileService> fs;

  TH1F* h_stripsPerNBX;
  TH1F* h_digiPerNBX;
  TH1F* h_padPerNBX;
  TH1F* h_coPadPerNBX;
  int m_bxnumber;
  long m_nevents;
  const int m_maxbx = 5;

  int m_digiPerNBX[2][2][36][10][384];
  int m_padPerNBX[2][2][36][10][384];
  int m_coPadPerNBX[2][2][36][10][384];

  TH1F* h_nstrips;
  TH1F* h_npads;
  TH1F* h_ncopads;
  
};
HitRateAnalysis::HitRateAnalysis(const edm::ParameterSet& iConfig)
{
  gemDigiInput_ = iConfig.getParameter<edm::InputTag>("gemDigiInput");
  gemPadDigiInput_ = iConfig.getParameter<edm::InputTag>("gemPadDigiInput");
  gemCoPadDigiInput_ = iConfig.getParameter<edm::InputTag>("gemCoPadDigiInput");

  m_bxnumber = 0;
  m_gemstrips.clear();
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int c = 0; c < 36; ++c)
	for (int f = 0; f < 10; ++f)
	  for (int k = 0; k < 384; ++k){
	    m_digiPerNBX[i][j][c][f][k] = 0;
	    m_padPerNBX[i][j][c][f][k] = 0;
	    m_coPadPerNBX[i][j][c][f][k] = 0;
	  }

  h_nstrips=fs->make<TH1F>("nstrips","nstrips",500,0,500);
  h_nstrips->GetYaxis()->SetTitle("Fraction");

  h_npads=fs->make<TH1F>("npads","npads",500,0,500);
  h_npads->GetYaxis()->SetTitle("Fraction");

  h_ncopads=fs->make<TH1F>("ncopads","ncopads",500,0,500);
  h_ncopads->GetYaxis()->SetTitle("Fraction");
  
  h_digiPerNBX=fs->make<TH1F>("digiPerNBX","digiPerNBX",10,0,10);
  h_digiPerNBX->GetXaxis()->SetTitle(Form("same digi hits within %2i events", m_maxbx));
  h_digiPerNBX->GetYaxis()->SetTitle("Fraction");

  h_stripsPerNBX=fs->make<TH1F>("stripsPerNBX","stripsPerNBX",10,0,10);
  h_stripsPerNBX->GetXaxis()->SetTitle(Form("same strips hits within %2i events", m_maxbx));
  h_stripsPerNBX->GetYaxis()->SetTitle("Fraction");
  
  h_padPerNBX=fs->make<TH1F>("padPerNBX","padPerNBX",10,0,10);
  h_padPerNBX->GetXaxis()->SetTitle(Form("same pad hits within %2i events", m_maxbx));
  h_padPerNBX->GetYaxis()->SetTitle("Fraction");
  
  h_coPadPerNBX=fs->make<TH1F>("coPadPerNBX","coPadPerNBX",10,0,10);
  h_coPadPerNBX->GetXaxis()->SetTitle(Form("same coPad hits within %2i events", m_maxbx));
  h_coPadPerNBX->GetYaxis()->SetTitle("Fraction");
  m_nevents = 0;
}
HitRateAnalysis::~HitRateAnalysis()
{
  cout << "m_nevents = " << m_nevents << endl;
  h_nstrips->Scale(1./m_nevents);
  h_npads->Scale(1./m_nevents);
  h_ncopads->Scale(1./m_nevents);
  h_stripsPerNBX->Scale(1./h_digiPerNBX->Integral());
  h_digiPerNBX->Scale(1./h_digiPerNBX->Integral());
  h_padPerNBX->Scale(1./h_padPerNBX->Integral());
  h_coPadPerNBX->Scale(1./h_coPadPerNBX->Integral());
}
void
HitRateAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GEMDigiCollection> gem_digis;  
  edm::Handle<GEMCSCPadDigiCollection> gemcscpad_digis;
  edm::Handle<GEMCSCCoPadDigiCollection> gemcsccopad_digis;
  iEvent.getByLabel(gemDigiInput_, gem_digis);
  iEvent.getByLabel(gemPadDigiInput_, gemcscpad_digis);
  iEvent.getByLabel(gemCoPadDigiInput_, gemcsccopad_digis);
  m_nevents++;
    
  for(GEMDigiCollection::DigiRangeIterator cItr = gem_digis->begin(); cItr != gem_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1){
	h_nstrips->Fill(digiIt->strip()-1);
	m_gemstrips.push_back(make_pair(id, digiIt->strip()-1));
	int region = 0;
	if (id.region() > 0) region = 1;
	m_digiPerNBX[region][id.layer()-1][id.chamber()-1][id.roll()-1][digiIt->strip()-1]++;
	if (36+id.region()*id.chamber() < 0) cout <<"error 1"<<endl;
	if (36+id.region()*id.chamber() >= 73) cout <<"error 2"<<endl;
	if (36+id.region()*id.chamber() == 36) cout <<"error 3"<<endl;
	if (id.layer()-1 < 0) cout <<"error 4"<<endl;
	if (id.layer()-1 >= 2) cout <<"error 5"<<endl;
	if (id.roll()-1 < 0) cout <<"error 6"<<endl;
	if (id.roll()-1 >= 10) cout <<"error 7"<<endl;
	if (digiIt->strip()-1 < 0) cout <<"error 8"<<endl;
	if (digiIt->strip()-1 >= 384){
	  cout <<"error 9 "<< digiIt->strip() <<endl;
	}
      }
    }
  }
  
  for(GEMCSCPadDigiCollection::DigiRangeIterator cItr = gemcscpad_digis->begin(); cItr != gemcscpad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMCSCPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1){
	h_npads->Fill(digiIt->pad()-1);
	int region = 0;
	if (id.region() > 0) region = 1;
	m_padPerNBX[region][id.layer()-1][id.chamber()-1][id.roll()-1][digiIt->pad()-1]++;
	if (digiIt->pad()-1 < 0) cout <<"error 10"<<endl;
	if (digiIt->pad()-1 >= 192){
	  cout <<"error 11 "<< digiIt->pad() <<endl;
	}	
      }
    }
  }

  for(GEMCSCCoPadDigiCollection::DigiRangeIterator cItr = gemcsccopad_digis->begin(); cItr != gemcsccopad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMCSCCoPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 1){
	h_ncopads->Fill(digiIt->first().pad()-1);
	int region = 0;
	if (id.region() > 0) region = 1;
	m_coPadPerNBX[region][id.layer()-1][id.chamber()-1][id.roll()-1][digiIt->first().pad()-1]++;
	if (digiIt->first().pad()-1 < 0) cout <<"error 12"<<endl;
	if (digiIt->first().pad()-1 >= 192){
	  cout <<"error 13 "<< digiIt->first().pad() <<endl;
	}
	//cout <<"id.layer() " << id.layer()-1 << endl;
      }
    }
  }
  m_bxnumber++;
  
  if (m_bxnumber == m_maxbx){
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
	for (int c = 0; c < 36; ++c)
	  for (int f = 0; f < 10; ++f)
	    for (int k = 0; k < 384; ++k){
	      int region = -1;
	      if (i == 1) region = 1;
	      GEMDetId id = GEMDetId(region, 1,1,j+1,c+1,f+1);
	      int nhits = 0;
	      for (unsigned int ns = 0; ns < m_gemstrips.size() ; ns++){
		if (m_gemstrips[ns].first == id){
		  if (m_gemstrips[ns].second == k){
		    nhits++;
		  }
		}
	      }
	      h_stripsPerNBX->Fill(nhits);

	      h_digiPerNBX->Fill(m_digiPerNBX[i][j][c][f][k]);

	      if (k < 192){
		h_padPerNBX->Fill(m_padPerNBX[i][j][c][f][k]);
		if (j == 0){
		  h_coPadPerNBX->Fill(m_coPadPerNBX[i][j][c][f][k]);
		}
	      }
	    }

    // reset all
    m_bxnumber = 0;
    m_gemstrips.clear();
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
	for (int c = 0; c < 36; ++c)
	  for (int f = 0; f < 10; ++f)
	    for (int k = 0; k < 384; ++k){
	      m_digiPerNBX[i][j][c][f][k] = 0;
	      m_padPerNBX[i][j][c][f][k] = 0;
	      m_coPadPerNBX[i][j][c][f][k] = 0;
	    }
    
  }
  
}

void HitRateAnalysis::beginJob(){}
void HitRateAnalysis::endJob(){}

//define this as a plug-in
DEFINE_FWK_MODULE(HitRateAnalysis);
