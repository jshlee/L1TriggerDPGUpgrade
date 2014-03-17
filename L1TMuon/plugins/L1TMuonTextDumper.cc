//////////////////////////////////////////////////////////////
// Upgraded Encdap Muon Track Finding Algorithm		    //
//							    //
// Info: A human-readable version of the firmware based     //
//       track finding algorithm which will be implemented  //
//       in the upgraded endcaps of CMS. DT and RPC inputs  //
//	 are not considered in this algorithm.      	    //
/// 							    //
// Author: M. Carver (UF)				    //
//////////////////////////////////////////////////////////////


#include "L1TriggerDPGUpgrade/L1TMuon/plugins/L1TMuonTextDumper.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h"
#include "L1Trigger/CSCTrackFinder/test/src/RefTrack.h"

#include <L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h>
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include <L1Trigger/CSCTrackFinder/interface/CSCTrackFinderDataTypes.h>

#include "DataFormats/MuonDetId/interface/GEMDetId.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PrimitiveConverter.h"
#include "BXAnalyzer.h"
#include "ZoneCreation.h"
#include "PatternRecognition.h"
#include "SortSector.h"
#include "Matching.h"
#include "Deltas.h"
#include "BestTracks.h"


using namespace L1TMuon;
using namespace std;

L1TMuonTextDumper::L1TMuonTextDumper(const PSet& p) {
  if( (_dogen = p.getUntrackedParameter<bool>("doGen",false)) ) {
    _geninput = p.getUntrackedParameter<edm::InputTag>("genSrc");
  }
  _tpinputs = p.getParameter<vector<edm::InputTag> >("primitiveSrcs");
  _convTrkInputs = 
    p.getParameter<vector<edm::InputTag> >("converterSrcs");
    
  LUTparam = p.getParameter<edm::ParameterSet>("lutParam");
    
  produces<L1TMuon::InternalTrackCollection> ("DataITC").setBranchAlias("DataITC");
}


void L1TMuonTextDumper::produce(edm::Event& ev, 
				const edm::EventSetup& es) {
			       
 		
  cout<<"Start TextDumper Producer::::: event = "<<ev.id().event()<<"\n\n";
  
  fprintf (write,"12345\n"); //<-- part of printing text file to send verilog code, not needed if George's package is included
  
  
  auto_ptr<L1TMuon::InternalTrackCollection> FoundTracks (new L1TMuon::InternalTrackCollection);
  
  vector<BTrack> PTracks[12];
 
  vector<TriggerPrimitiveRef> tester;
  //vector<InternalTrack> FoundTracks;
  
  //////////////////////////////////////////////
  ////////// Get Generated Muons ///////////////
  //////////////////////////////////////////////
  
  edm::Handle<vector<reco::GenParticle>> GenMuons;
  vector<reco::GenParticle>::const_iterator GI;
  ev.getByLabel("genParticles",GenMuons);
  bool gpir = false, endcap1 = false, endcap2 = false;
  int etaindex = -99;
  reco::GenParticle GeneratorMuon;
  for(GI=GenMuons->begin();GI!=GenMuons->end();GI++){
  	
    const reco::GenParticle GenMuon = *GI;
    GeneratorMuon = GenMuon;
    double pt = GenMuon.pt(), eta = GenMuon.eta(), phi = GenMuon.phi(), mass = GenMuon.mass();
    int charge = GenMuon.charge();
	
    cout<<"Gen Particle Info::::\nPt = "<<pt<<", phi = "<<phi<<", eta = "<<eta<<", mass = "<<mass<<" and charge = "<<charge<<"\n\n";
  	
    if((fabs(eta) > 1.2) && (fabs(eta) <= 2.4) && (pt >= 5))
      gpir = true;
	
    if(eta > 0)
      endcap1 = true;
	
    if(eta < 0)
      endcap2 = true;
	
    for(int y=0;y<24;y++){
	
      double low = -2.4;
      int mult = y;
      if(y > 11){
	low = 1.2;
	mult = y - 12;
      }
			
	
      double el = low + 0.1*mult;
      double eh = low + 0.1*(mult + 1);
		
      if(eta >= el && eta < eh)
	etaindex = y;
	
    }
		
  }
  
  
  if(etaindex != -99)
    cout<<"\neta index = "<<etaindex<<"\n";//
  
  //////////////////////////////////////////////
  ///////// Get Trigger Primitives /////////////  Retrieve TriggerPrimitives from the event record
  //////////////////////////////////////////////
  
  auto tpsrc = _tpinputs.cbegin();
  auto tpend = _tpinputs.cend();
  for( ; tpsrc != tpend; ++tpsrc ) {
    edm::Handle<TriggerPrimitiveCollection> tps;
    ev.getByLabel(*tpsrc,tps);
    auto tp = tps->cbegin();
    auto tpend = tps->cend();

    for( ; tp != tpend; ++tp ) {
      std::cout << "jason: tp->subsystem()" << tp->subsystem() << std::endl;
      TriggerPrimitiveRef tpref(tps,tp - tps -> cbegin());
		
      tester.push_back(tpref);
		
      if(tp->subsystem() == 1){
	cout<<"\ntrigger prim found station:"<<tp->detId<CSCDetId>().station()<<endl;
		
	if((tp->detId<CSCDetId>().station() == 4) && (fabs(GeneratorMuon.eta()) < 1.7)){
		
	  if(tp->detId<CSCDetId>().endcap() == 1)
	    ME42test1->Fill(GeneratorMuon.phi());
			
	  if(tp->detId<CSCDetId>().endcap() == 2)
	    ME42test2->Fill(GeneratorMuon.phi());
		
	}
	  
	if((tp->detId<CSCDetId>().station() == 1) && (tp->detId<CSCDetId>().ring() == 4))
	  ME1gangnedtest->Fill(tp->getCSCData().strip);
			
	if((tp->detId<CSCDetId>().station() == 1) && (tp->detId<CSCDetId>().ring() == 1))
	  ME11gangnedtest->Fill(tp->getCSCData().strip);
	  
	if(tp->detId<CSCDetId>().triggerSector() == 1){
		
	  cout<<"there are sector 1 csc hits\n\n";
	  cout<<"strip = "<<tp->getCSCData().strip<<endl;
	  striph->Fill(tp->getCSCData().strip);
	}
      }
      if(tp->subsystem() == 3){
	if((tp->detId<GEMDetId>().station() == 1))
	  h_GE11->Fill(1);			
      }
    }
  }
  vector<ConvertedHit> CHits[12];
  MatchingOutput MO[12];
 
  for(int SectIndex=0;SectIndex<12;SectIndex++){//perform TF on all 12 sectors

 
 
    //////////////////////////////////////////////////////  Input is raw hit information from 
    ///////////////// TP Conversion //////////////////////  Output is vector of Converted Hits
    //////////////////////////////////////////////////////


    vector<ConvertedHit> ConvHits = PrimConv(tester,SectIndex);
    CHits[SectIndex] = ConvHits;
	
    for(vector<ConvertedHit>::iterator i=ConvHits.begin();i!=ConvHits.end();i++){
	
      if(i->TP()->detId<CSCDetId>().ring() == 4)
	ME1gangnedtest->Fill(i->Phi());
	
    }
   
 
  
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
    ////////////////////////////////print values for input into Alex's emulator code/////////////////////////////////////////////////////
    //for(vector<ConvertedHit>::iterator h = ConvHits.begin();h != ConvHits.end();h++){
	
    //if((h->Id()) > 9){h->SetId(h->Id() - 9);h->SetStrip(h->Strip() + 128);}	
    //fprintf (write,"0	1	1 	%d	%d\n",h->Sub(),h->Station());
    //fprintf (write,"1	%d	%d 	%d\n",h->Quality(),h->Pattern(),h->Wire());
    //fprintf (write,"%d	0	%d\n",h->Id(),h->Strip());	
    //}
    ////////////////////////////////print values for input into Alex's emulator code/////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  


    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////  Takes the vector of converted hits and groups into 3 groups of hits
    ////////////////////// BX Grouper ////////////////////  which are 3 BX's wide. Effectively looking 2 BX's into the future and 
    //////////////////////////////////////////////////////  past from the central BX, this analyzes a total of 5 BX's.
    //////////////////////////////////////////////////////
 
 
    vector<vector<ConvertedHit>> GroupedHits = GroupBX(ConvHits);
 
 
    ////////////////////////////////////////////////////////  Creates a zone for each of the three groups created in the BX Grouper module.
    ////////// Creat Zones for pattern Recognition /////////  The output of this module not only contains the zones but also the 
    ////////////////////////////////////////////////////////  reference back to the TriggerPrimitives that went into making them.
	
    vector<ZonesOutput> Zout = Zones(GroupedHits);
   

    ///////////////////////////////
    ///// Pattern Recognition /////  Applies pattern recognition logic on each of the 3 BX groups and assigns a quality to each keystrip in the zone.
    ///// & quality assinment /////  The delete duplicate patterns function looks at the 3 BX groups and deletes duplicate patterns found from the
    ///////////////////////////////  same hits. This is where the BX analysis ends; Only 1 list of found patterns is given to the next module.
  

    vector<PatternOutput> Pout = Patterns(Zout);
  
    PatternOutput Test = DeleteDuplicatePatterns(Pout);
 
    cout << "test PatternOutput " << endl;
    PrintQuality(Test.detected);
 

    ///////////////////////////////
    //////Sector Sorting/////////// Sorts through the patterns found in each zone and selects the best three per zone to send to the next module.
    ///////Finding 3 Best Pattern// 
    ///////////////////////////////
  
  
    SortingOutput Sout = SortSect(Test);
	
 
    //////////////////////////////////
    ///////// Match ph patterns ////// Loops over each sorted pattern and then loops over all possible triggerprimitives which could have made the pattern
    ////// to segment inputs ///////// and matches the associated full precision triggerprimitives to the detected pattern. 
    //////////////////////////////////   
      

    MatchingOutput Mout = PhiMatching(Sout);
    MO[SectIndex] = Mout;

    /////////////////////////////////
    //////// Calculate delta //////// Once we have matched the hits we calculate the delta phi and theta between all 
    ////////    ph and th    //////// stations present. 
    /////////////////////////////////
  
  
    vector<vector<DeltaOutput>> Dout = CalcDeltas(Mout);////
 

    /////////////////////////////////
    /////// Sorts and gives /////////  Loops over all of the found tracks(looking across zones) and selects the best three per sector. 
    ////// Best 3 tracks/sector /////  Here ghost busting is done to delete tracks which are comprised of the same associated stubs. 
    /////////////////////////////////  
  
  
    vector<BTrack> Bout = BestTracks(Dout);
    PTracks[SectIndex] = Bout;
   
  
  	
  }
 
  ////////////////////////////////////
  //// Ghost Cancellation between ////not done correct
  //////////   sectors     ///////////
  ////////////////////////////////////
 
  for(int i1=0;i1<36;i1++){
 
    for(int i2=i1+1;i2<36;i2++){
		
      int sec[2] = {i1/3,i2/3};
      int num[2] = {i1%3,i2%3};
		
      bool wrap = ((sec[0] == 0 || sec[0] == 6) && (fabs(sec[0] - sec[1]) == 5));
      bool same = (sec[0] == sec[1]);
      bool toofar = (fabs(sec[0] - sec[1]) > 1);

      if((same || toofar) && !wrap)//if same chamber or more than one chamber away dont do. !wrap allows comparison between sectors 0&5 and 6&11
	continue;
			
      if(!PTracks[sec[0]][num[0]].AHits.size() || !PTracks[sec[1]][num[1]].AHits.size())
	continue;
		
      int sh_seg = 0;
		
      cout<<"\nComparing adjacent sectors\n";

      for(int sta=0;sta<4;sta++){//the part which is done incorrectly
		
	if((PTracks[sec[0]][num[0]].AHits[sta].Phi() == -999) || (PTracks[sec[1]][num[1]].AHits[sta].Phi() == -999))
	  continue;
				
		
	if((PTracks[sec[0]][num[0]].AHits[sta].Id() == PTracks[sec[1]][num[1]].AHits[sta].Id())
	   && (PTracks[sec[0]][num[0]].AHits[sta].Strip() == PTracks[sec[1]][num[1]].AHits[sta].Strip())
	   && (PTracks[sec[0]][num[0]].AHits[sta].Wire() == PTracks[sec[1]][num[1]].AHits[sta].Wire())){
				
	  sh_seg++;
	}
			
      }
		
      if(sh_seg){//if any segments are shared delete the track with lower rank
			
	BTrack tmp;//default null track to replace ghost with
			
	if(PTracks[sec[0]][num[0]].winner.Rank() >= PTracks[sec[1]][num[1]].winner.Rank())
	  PTracks[sec[1]][num[1]] = tmp;
	else
	  PTracks[sec[0]][num[0]] = tmp;
		
      }
    }
  }
 
 
  ////////////////////////////////////
  /// Sorting through all sectors ////
  ///   to find 4 best muons      ////
  ////////////////////////////////////
 
 
  BTrack FourBest[4];
  vector<BTrack> PTemp[12] = PTracks;
  int windex[4] = {-1,-1,-1,-1};
 
 
 
  for(int i=0;i<4;i++){
 
    for(int j=0;j<36;j++){
 
		
      if(!PTemp[j/3][j%3].phi)//no track
	continue;
			
      if((windex[0] == j) || (windex[1] == j) || (windex[2] == j) || (windex[3] == j))//already picked
	continue;
		
      if(PTracks[j/3][j%3].winner.Rank() > FourBest[i].winner.Rank()){
			
	FourBest[i] = PTemp[j/3][j%3];
	windex[i] = j;
				
      } 

    }
  }
  ///////////////////////////////////
  /// Make Internal track if ////////
  /////// tracks are found //////////
  ///////////////////////////////////

  bool epir = false;
  
  for(int fbest=0;fbest<4;fbest++){
  
    if(FourBest[fbest].phi){
	
      epir = true;
	
      InternalTrack tempTrack;
      tempTrack.setType(2); 
	
      tempTrack.phi = FourBest[fbest].phi;
      tempTrack.theta = FourBest[fbest].theta;
      tempTrack.rank = FourBest[fbest].winner.Rank();
      tempTrack.deltas = FourBest[fbest].deltas;
		
      for(vector<ConvertedHit>::iterator A = FourBest[fbest].AHits.begin();A != FourBest[fbest].AHits.end();A++){
		
	if(A->Phi() != -999){
			
	  tempTrack.addStub(A->TP());
	  //cout<<"Q: "<<A->Quality()<<", keywire: "<<A->Wire()<<", strip: "<<A->Strip()<<endl;
	}
			
      }
		
      FoundTracks->push_back(tempTrack);
    }
  }
  
  ///////////////////////////////////////////////
  //// Pt assignment //
  ///////////////////////////////////////////////
  edm::ESHandle< L1MuTriggerScales > scales ;
  es.get< L1MuTriggerScalesRcd >().get( scales ) ;

  edm::ESHandle< L1MuTriggerPtScale > ptScale ;
  es.get< L1MuTriggerPtScaleRcd >().get( ptScale ) ;

  //  edm::ParameterSet ptLUTset = p.getParameter<edm::ParameterSet>("PTLUT");
  CSCTFPtLUT* ptLUT_ = new CSCTFPtLUT(LUTparam, scales.product(),ptScale.product());
   
  cout << "ptLUT_ " << ptLUT_ << endl;
  cout << "FoundTracks->size() " << FoundTracks->size() << endl;
  //  for (auto foundTrack : FoundTracks){

  auto stub = FoundTracks->begin();
  auto stend = FoundTracks->end();
  for( ; stub != stend; ++stub ) {

    // cout << "foundTrack->phiValue() " << stub->phiValue() << endl;
    // cout << "foundTrack->etaValue() " << stub->etaValue() << endl;
    // cout << "foundTrack->ptValue() " << stub->ptValue() << endl;

    // ptadd address;
    // address.delta_phi_12 = delta_phi_12;
    // address.delta_phi_23 = delta_phi_23;
    // address.track_eta = track_eta;
    // address.track_mode = track_mode;
    // address.track_fr = track_fr;
    // address.delta_phi_sign = delta_phi_sign;
    
    // vector<csc::L1Track>::iterator titr = tftks.begin();

    // for(; titr != tftks.end(); titr++){
    //   ptadd thePtAddress(titr->ptLUTAddress());
    //   ptdat thePtData = ptLUT_->Pt(thePtAddress);

    //   if (thePtAddress.track_fr){
    //     titr->setRank(thePtData.front_rank);
    //     titr->setChargeValidPacked(thePtData.charge_valid_front);
    //   }
    //   else {
    //     titr->setRank(thePtData.rear_rank);
    //     titr->setChargeValidPacked(thePtData.charge_valid_rear);
    //   }
    //   if ( ((titr->ptLUTAddress()>>16)&0xf)==15 ){
    //     int unmodBx = titr->bx();
    //     titr->setBx(unmodBx+2);
    //   }
    // }
  }
  ///////////////////////////////////////////////
  //// Below here is working additions to make //
  //// efficiency plots and can be neglected ////
  ///////////////////////////////////////////////
  
  if(gpir)
    eff->Fill(FoundTracks->size());
  if(gpir && endcap2)
    eff2->Fill(FoundTracks->size());

  if(gpir && epir && etaindex != -99)
    fpire[etaindex]++;

  bool detectorinefficiency = false;
  
  if(gpir && !epir){
  
    if(endcap1)
      trigprimsize->Fill(tester.size());
    else if(endcap2)
      trigprimsize2->Fill(tester.size());
	
    if(tester.size() < 2)
      detectorinefficiency = true;
	
    if(tester.size() >= 2){
	
      int contribution = 0;
      for(vector<TriggerPrimitiveRef>::iterator C1 = tester.begin();C1 != tester.end();C1++){
	if ((*C1)->subsystem() == TriggerPrimitive::kCSC){
	  int station = (*C1)->detId<CSCDetId>().station();
	  switch(station){
	  case(1):contribution |= 8;break;
	  case(2):contribution |= 4;break;
	  case(3):contribution |= 2;break;
	  case(4):contribution |= 1;break;
	  default:cout<<"Station is out of range\n";
	  }
	}
      }
		
      if(endcap1)
	st_cont->Fill(contribution);
      else if(endcap2)
	st_cont2->Fill(contribution);
		
      if(contribution == 1 || contribution == 2 || contribution == 3 || contribution == 4 || contribution == 8)
	detectorinefficiency = true;
		
      if(contribution == 6){//had stations 2&3 but no track
	cout<<"Stations 2 & 3\nevent = "<<ev.id().event()<<"\nSectIndex = "<<
	  ((tester[0]->detId<CSCDetId>().endcap()) - 1)*6 + (tester[0]->detId<CSCDetId>().triggerSector()) - 1<<"\n\n";
      }
		
      if(contribution == 11){//had stations 1&3&4 but no track
	cout<<"Stations 1 & 3 & 4\nevent = "<<ev.id().event()<<"\nSectIndex = "<<
	  ((tester[0]->detId<CSCDetId>().endcap()) - 1)*6 + (tester[0]->detId<CSCDetId>().triggerSector()) - 1<<"\n\n";
      }
		
      if(contribution == 10){//had stations 1&3 but no track
	cout<<"Stations 1 & 3\nevent = "<<ev.id().event()<<"\nSectIndex = "<<
	  ((tester[0]->detId<CSCDetId>().endcap()) - 1)*6 + (tester[0]->detId<CSCDetId>().triggerSector()) - 1<<"\n\n";
      }
		
		
		
      if(contribution == 12){//had stations 1 & 2 but didn't make a track
		
	int sindex0 = ((tester[0]->detId<CSCDetId>().endcap()) - 1)*6 + (tester[0]->detId<CSCDetId>().triggerSector()) - 1;
	int sindex1 = ((tester[1]->detId<CSCDetId>().endcap()) - 1)*6 + (tester[1]->detId<CSCDetId>().triggerSector()) - 1;
	int zhit1 = CHits[sindex0][0].Zhit(), zhit2 = CHits[sindex0][1].Zhit();
			
	bool ptf = (fabs(zhit2 - zhit1) > 15);
	bool thwindow = false;
			
	for(int zone=0;zone<4;zone++){
	  for(int winner=0;winner<3;winner++){
			
	    ThOutput thmatch = MO[sindex0].ThetaMatch();
	    int dth[6][4] = {{-999,-999,-999,-999},{-999,-999,-999,-999},{-999,-999,-999,-999},{-999,-999,-999,-999},{-999,-999,-999,-999},{-999,-999,-999,-999}};
	    for(int s1=0;s1<3;s1++){
	
	      for(int s2=s1+1;s2<4;s2++){
			
			
		///////////////////////// There is a further index on dTh because there are 4 dth combinations 
		/// calc delta theta  /// possible if there are two theta segments for both stations. 
		///////////////////////// EXPLAIN ABOUT [I+J] AND [I+J+1] 
					
		for(int i=0;i<2;i++){
						
					
		  for(int j=0;j<2;j++){
						
		    int thi = thmatch[zone][winner][s1][i].Theta();
		    int thj = thmatch[zone][winner][s2][j].Theta();
		    int deltath = thi - thj;
					
					
		    if((s1 == 0) && (thi != -999) && (thj != -999)){///need to fix still////
								
		      if(!i){dth[s2-1][i+j] = deltath;}
		      else{dth[s2-1][i+j+1] = deltath;}
								
								
		    }
		    else if((s1 != 0) && (thi != -999) && (thj != -999)){
							
		      if(!i){dth[s1+s2][i+j] = deltath;}
		      else{dth[s1+s2][i+j+1] = deltath;}
		    }
							
		  }
		}
		
	      }
	    }
			
	    for(int b1=0;b1<6;b1++){
	      for(int b2=0;b2<2;b2++){
					
		if((dth[b1][b2] != -999) && (fabs(dth[b1][b2]) <= 4))
		  thwindow = true;
	      }
	    }
	  }
	}
			
			
	if(sindex0 != sindex1){
			
	  if((sindex0 > 5) || (sindex1 > 5))
	    st12errors->Fill(1);
	  else
	    st12errors->Fill(0);
	}
	else if(sindex0 == sindex1){
	  bool samezone = false;
			
	  for(unsigned int i1 = 0;i1 != tester.size();i1++){
			
	    for(unsigned int i2 = i1+1;i2 != tester.size();i2++){
				
					
	      int s1 = tester[i1]->detId<CSCDetId>().station(), s2 = tester[i2]->detId<CSCDetId>().station();
	      vector<int> z1 = CHits[sindex0][i1].ZoneContribution(), z2 = CHits[sindex0][i2].ZoneContribution();
					
	      if(s1 != s2){
					
		for(vector<int>::iterator a1=z1.begin();a1 != z1.end();a1++){
		  for(vector<int>::iterator a2=z2.begin();a2 != z2.end();a2++){
							
		    if((*a1) == (*a2))
		      samezone = true;
								
		  }
		}
	      }
	    }
	  }
			
			
			
	  if(!samezone){
				
	    if(sindex0 > 5)
	      st12errors->Fill(3);
	    else
	      st12errors->Fill(2);
					
	  }
	  else if(samezone){
				
					
	    if(ptf){
						
	      if(sindex0 > 5)
		st12errors->Fill(5);
	      else
		st12errors->Fill(4);
	    }
	    else if(!ptf){
						
	      if(!thwindow){
						
		if(sindex0 > 5)
		  st12errors->Fill(7);
		else
		  st12errors->Fill(6);
								
		cout<<"\n!thwindow\n";
	      }
	      //else if(thwindow){
						
						
							
						
	      //}
	    }
	  }
	}
      }
		
      if((contribution > 4) && (contribution != 8)){
		
	cout<<"\nMISSED\n";
			
	if(contribution == 12)
	  cout<<"ONLY STATIONS 1 AND 2 : (\nevent = "<<ev.id().event()<<"\n\n";
		
	MissVsEta->Fill(GeneratorMuon.eta());
	MissVsPhi->Fill(GeneratorMuon.phi());
	MissVsPt->Fill(GeneratorMuon.pt());	
		
      }
    }
  }
  
  
  if(gpir && !detectorinefficiency)
    gpire[etaindex]++;
  

  if(gpir && (FoundTracks->size() == 1)){
  
    cout<<"\nFOUND ONE MUON-------------Sector "<<windex[0]/3<<"\n";
	
    int ecap = 0;
	
    vector<ConvertedHit> ahits = FourBest[0].AHits;
	
    cout<<"ahits.size() = "<<ahits.size()<<endl;
    cout<<"ahits[0].Phi() "<< ahits[0].Phi() << endl;
    cout<<"ahits[1].Phi() "<< ahits[1].Phi() << endl;
    cout<<"ahits[2].Phi() "<< ahits[2].Phi() << endl;
    cout<<"ahits[3].Phi() "<< ahits[3].Phi() << endl;
	
    
    if(FourBest[0].AHits[0].Phi() != -999){
      if(FourBest[0].AHits[0].TP()->detId<CSCDetId>().endcap()){
	ecap = FourBest[0].AHits[0].TP()->detId<CSCDetId>().endcap();
	cout<<"\n1\n";}
    }
    else if(FourBest[0].AHits[1].Phi() != -999){
      if(FourBest[0].AHits[1].TP()->detId<CSCDetId>().endcap()){
	ecap = FourBest[0].AHits[1].TP()->detId<CSCDetId>().endcap();
	cout<<"\n2\n";}
    }
    else if(FourBest[0].AHits[2].Phi() != -999){
      if(FourBest[0].AHits[2].TP()->detId<CSCDetId>().endcap()){
	ecap = FourBest[0].AHits[2].TP()->detId<CSCDetId>().endcap();
	cout<<"\n3\n";}
    }
    else if(FourBest[0].AHits[3].Phi() != -999){
      if(FourBest[0].AHits[3].TP()->detId<CSCDetId>().endcap()){
	ecap = FourBest[0].AHits[3].TP()->detId<CSCDetId>().endcap();
	cout<<"\n4\n";}
    }
		
    cout<<"ecap "<< ecap << endl;

    cout<<"\n5\n";
	
    int cont = 0, numTP = 0;
    for(vector<TriggerPrimitiveRef>::iterator C1 = tester.begin();C1 != tester.end();C1++){
			
      int stat = 0;
      cout<<"\n2\n";
      
      if ((*C1)->subsystem() == TriggerPrimitive::kCSC){
	if((*C1)->detId<CSCDetId>().endcap() != ecap){
	  stat = (*C1)->detId<CSCDetId>().station();
	  numTP++;
	}
      }
      cout<<"\n3\n";
		
      switch(stat){
			
	//case(0):cont |= 0;break
      case(1):cont |= 8;break;
      case(2):cont |= 4;break;
      case(3):cont |= 2;break;
      case(4):cont |= 1;break;
      default:cout<<"Station is out of range\n";
			
      }
    }
	
    if(numTP < 2){
      detectorineff->Fill(numTP);
    }
    else{
	
      switch(cont){
      case(1):detectorineff->Fill(2);break;//only st 4 present
      case(2):detectorineff->Fill(3);break;//only st 3 present
      case(3):detectorineff->Fill(4);break;//only st's 3 and 4 present
      case(4):detectorineff->Fill(5);break;//only st 2 present
      case(8):detectorineff->Fill(6);break;//only st 1 present
      }	
    }	
  
  }
  
  if(gpir && (FoundTracks->size() == 2)){
  
  	
  
    int sectors[2] = {windex[0]/3,windex[1]/3};
    sector1->Fill(sectors[0]);
    sector2->Fill(sectors[1]);
    secdiff->Fill(fabs(sectors[0] - sectors[1]));
	
    if(fabs(sectors[0] - sectors[1]) == 0)
      cout<<"\nTWO IN SAME SECTOR\n";
	
    cout<<"\nTWO MUONS FOUND------- Sectors: "<<sectors[0]<<" and "<<sectors[1]<<"\n";
	
  
  }
  
  if(gpir && (FoundTracks->size() == 3)){
  
    cout<<"\nFOUND THREE MUONS---------- Sectors: "<<windex[0]/3<<", "<<windex[1]/3<<" and "<<windex[2]/3<<"\n";
  }

  
  //  cout<<"Begin Put function\n\n";
  ev.put( FoundTracks, "DataITC");
  cout<<"End TextDump Prducer:::::::::::::::::::::::::::\n:::::::::::::::::::::::::::::::::::::::::::::::::\n\n";

}//analyzer

void L1TMuonTextDumper::beginJob()
{

  //cout<<"Begin TextDump Prducer:::::::::::::::::::::::::::\n:::::::::::::::::::::::::::::::::::::::::::::::::\n\n";
  ///////////////////////////
  ////// Histogram //////////
  ////// Declaration ////////
  ///////////////////////////
	
  TFileDirectory dir = histofile->mkdir("1");//
  TFileDirectory dir1 = dir.mkdir("2");
  TFileDirectory dir2 = dir.mkdir("3");
	
  ///////////////////////////
  /////// Output ////////////
  ///// Text Files //////////
  ///////////////////////////
	
	
  write = fopen ("zone0.txt","w");
  dphi = fopen("dphi1.txt","w");
  tptest = fopen("dth.txt","w");
	
  h_GE11 = dir1.make<TH1F>("GE11","GE11",10,0,10);h_GE11->SetFillColor(2);
  striph = dir1.make<TH1F>("striph","TP strip distribution",250,0,250);striph->SetFillColor(2);
  eff = dir.make<TH1F>("eff","If GenParticle how many EmulatorMuons (Endcap 1)",5,0,5);eff->SetFillColor(4);
  eff2 = dir.make<TH1F>("eff2","If GenParticle how many EmulatorMuons (Endcap 2)",5,0,5);eff2->SetFillColor(4);
  trigprimsize = dir.make<TH1F>("trigprimsize","How many TP's in event if no found Muon (Endcap 1)",9,0,9);trigprimsize->SetFillColor(4);
  trigprimsize2 = dir.make<TH1F>("trigprimsize2","How many TP's in event if no found Muon (Endcap 2)",9,0,9);trigprimsize2->SetFillColor(4);
  st_cont = dir.make<TH1F>("st_cont","Stations Present if no found Muon and 2 or more TP's (Endcap 1)",16,0,16);st_cont->SetFillColor(4);
  st_cont2 = dir.make<TH1F>("st_cont2","Stations Present if no found Muon and 2 or more TP's (Endcap 2)",16,0,16);st_cont2->SetFillColor(4);
  sector1 = dir.make<TH1F>("sector1","Sector 1 if 2 Muons found", 12,0,12);sector1->SetFillColor(4);
  sector2 = dir.make<TH1F>("sector2","Sector 2 if 2 Muons found", 12,0,12);sector2->SetFillColor(4);
  secdiff = dir.make<TH1F>("secdiff","Sector1 - Sector2 if 2 Muons found",12,0,12);secdiff->SetFillColor(4);
  MissVsEta = dir.make<TH1F>("MissVsEta","Eta Distribution of Missed Muons",50,-2.5,2.5);MissVsEta->SetFillColor(4);
  MissVsPhi = dir.make<TH1F>("MissVsPhi","Phi Distribution of Missed Muons",70,-3.5,3.5);MissVsPhi->SetFillColor(4);
  MissVsPt = dir.make<TH1F>("MissVsPt","Pt Distribution of Missed Muons",200,0,200);MissVsPt->SetFillColor(4);
	
  ME42test1 = dir.make<TH1F>("ME42test1","Phi Dist. of ME42 Hits (Endcap 1)",70,-3.5,3.5);ME42test1->SetFillColor(4);
  ME42test2 = dir.make<TH1F>("ME42test2","Phi Dist. of ME42 Hits (Endcap 2)",70,-3.5,3.5);ME42test2->SetFillColor(4);
  ME1gangnedtest = dir.make<TH1F>("ME1gangnedtest","Strip Dist. of ME11A hits",5000,0,5000);ME1gangnedtest->SetFillColor(4);
  ME11gangnedtest = dir.make<TH1F>("ME11gangnedtest","Strip Dist. of ME11 hits",200,0,200);ME11gangnedtest->SetFillColor(4);
	
  st12errors = dir.make<TH1F>("st12errors","Error Indication for a miss with TP's in station 1 and 2",8,0,8);st12errors->SetFillColor(4);
	
  detectorineff = dir.make<TH1F>("detectorineff","Detector Inefficiencies", 7,0,7);detectorineff->SetFillColor(4);
	
  for(int i=0;i<24;i++){
    gpire[i] = 0;
    fpire[i] = 0;
  }
	
	
	
}
void L1TMuonTextDumper::endJob()
{

  fclose (write);
  fclose (dphi);
  fclose (tptest);
  TFileDirectory dir = histofile->mkdir("1");
	
  cout<<"\n\n\nfpire = ";
  for(int i=0;i<24;i++)
    cout<<fpire[i]<<"   ";
	
  cout<<"\ngpire = ";
  for(int ii=0;ii<24;ii++)
    cout<<gpire[ii]<<"   ";
		
  cout<<"\nTHE END"<<endl;
	
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonTextDumper);
