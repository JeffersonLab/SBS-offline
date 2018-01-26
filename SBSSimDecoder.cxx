//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   SBSSimDecoder
//
//   Decoder for SBS simulation data
//
//   Interprets event buffer from input as SBSSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "SBSSimDecoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
//#include "TSBSDBManager.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
using namespace Podd;

//EFuchey: 2016/12/10: it is necessary to declare the SBSDBManager as a static instance here 
// (and not make it a member) because it is used by functions whic are defined as "static inline".
//static SBSDBManager* fManager = SBSDBManager::GetInstance();
static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in TreeSearch
enum EProjType { kUPlane = 0, kVPlane =1, kXPlane = 2, kYPlane = 3};

typedef vector<int>::size_type vsiz_t;

//-----------------------------------------------------------------------------
SBSSimDecoder::SBSSimDecoder()
{
  // Constructor

  fMCCherHits = new TClonesArray( "SBSSimPMTHit",      2000 );
  fMCCherClus = new TClonesArray( "SBSSimCherCluster", 500 );
  
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
  gSystem->Load("libsbscer.so");  // for library objects
  // maybe we can load dynamic libraries 
}

//-----------------------------------------------------------------------------
SBSSimDecoder::~SBSSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );
}

//-----------------------------------------------------------------------------
Int_t SBSSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define global variables for the MC quantities. Extends the base
  // class method.

  const char* const here = "SBSSimDecoder::DefineVariables";
  
  if( mode == THaAnalysisObject::kDefine && fIsSetup )
    return THaAnalysisObject::kOK;
  
  SimDecoder::DefineVariables( mode );
  
  
  //Cherenkov hits
  RVarDef vars[] = {
    // PMT hits
    { "ch.pmt.id",      "MC Cher PMT ID",             "fMCCherHits.SBSSimPMTHits.fID"          },
    { "ch.pmt.src",     "MC Cher PMT Source",         "fMCCherHits.SBSSimPMTHits.fSource"      },
    { "ch.pmt.type",    "MC Cher PMT Type",           "fMCCherHits.SBSSimPMTHits.fType"        },
    { "ch.pmt.MCID",    "MC Cher PMT MC track ID",    "fMCCherHits.SBSSimPMTHits.fMCtrackID"   },
    { "ch.pmt.MCPID",   "MC Cher PMT MC track PID",   "fMCCherHits.SBSSimPMTHits.fMCtrackPID"  },
    { "ch.pmt.origVol", "MC Cher PMT origin volume",  "fMCCherHits.SBSSimPMTHits.fOrigVolFlag" },
    { "ch.pmt.x",       "MC Cher PMT X (m)",          "fMCCherHits.SBSSimPMTHits.fXPMT"        },
    { "ch.pmt.y",       "MC Cher PMT Y (m)",          "fMCCherHits.SBSSimPMTHits.fYPMT"        },
    { "ch.pmt.Npe",     "MC Cher PMT N pe",           "fMCCherHits.SBSSimPMTHits.fNpe"         },
    { "ch.pmt.t0",      "MC Cher PMT t0 (ns)",        "fMCCherHits.SBSSimPMTHits.fTime"        },
    { "ch.pmt.tr",      "MC Cher PMT rise time (ns)", "fMCCherHits.SBSSimPMTHits.fTDCtime1"    },
    { "ch.pmt.tf",      "MC Cher PMT fall time (ns)", "fMCCherHits.SBSSimPMTHits.fTDCtime2"    },
    { "ch.pmt.detID",   "MC Cher PMT detector ID",    "fMCCherHits.SBSSimPMTHits.fDetID"       },
    { "ch.pmt.chan",    "MC Cher PMT chan",           "fMCCherHits.SBSSimPMTHits.fChannel"     },
    { "ch.pmt.row",     "MC Cher PMT row",            "fMCCherHits.SBSSimPMTHits.fPMTrow"      },
    { "ch.pmt.col",     "MC Cher PMT col",            "fMCCherHits.SBSSimPMTHits.fPMTcol"      },
    
    // PMT clusters
    { "ch.cl.size",    "MC Cher clus size",                "fMCCherClus.SBSSimCherCluster.fSize"           },
    { "ch.cl.X",       "MC Cher clus avg X (m)",           "fMCCherClus.SBSSimCherCluster.fXcenter"        },
    { "ch.cl.Y",       "MC Cher clus avg Y (m)",           "fMCCherClus.SBSSimCherCluster.fXcenter"        },
    { "ch.cl.X_w",     "MC Cher clus avg X (w Npe) (m)",   "fMCCherClus.SBSSimCherCluster.fXcenter_w"      },
    { "ch.cl.Y_w",     "MC Cher clus avg X (w Npe) (m)",   "fMCCherClus.SBSSimCherCluster.fYcenter_w"      },
    { "ch.cl.Npe",     "MC Cher clus Npe total",           "fMCCherClus.SBSSimCherCluster.fYcenter_w"      },
    { "ch.cl.tr_avg",  "MC Cher clus avg rise time (ns)",  "fMCCherClus.SBSSimCherCluster.fMeanRisingTime" },
    { "ch.cl.tr_rms",  "MC Cher clus rise time rms (ns)",  "fMCCherClus.SBSSimCherCluster.fRisingTimeRMS"  },
    { "ch.cl.tf_avg",  "MC Cher clus avg fall time (ns)",  "fMCCherClus.SBSSimCherCluster.fMeanRisingTime" },
    { "ch.cl.tf_rms",  "MC Cher clus fall time rms (ns)",  "fMCCherClus.SBSSimCherCluster.fRisingTimeRMS"  },
    { "ch.cl.MCID",    "MC Cher clus track G4ID",          "fMCCherClus.SBSSimCherCluster.fMCtrackID"      },
    { "ch.cl.MCPID",   "MC Cher clus track G4PID",         "fMCCherClus.SBSSimCherCluster.fMCtrackPID"     },
    /*
      if something is to be added, can be found easily in libsbsgem::SBSSimDecoder
    */
    { 0 }
  };

  return THaAnalysisObject::
    DefineVarsFromList( vars, THaAnalysisObject::kRVarDef,
			mode, "", this, MC_PREFIX, here );
}

//-----------------------------------------------------------------------------
void SBSSimDecoder::Clear( Option_t* opt )
{
  // Clear track and plane data

  SimDecoder::Clear(opt);   // clears fMCHits, fMCTracks and fMCPoints

  fPMTMap.clear();
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
int SBSSimDecoder::LoadEvent(const UInt_t* evbuffer )
#else
int SBSSimDecoder::LoadEvent(const Int_t* evbuffer )
#endif
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  int ret = DoLoadEvent( evbuffer );

  if( fDoBench ) fBench->Stop("physics_decode");

  return ret;
}

//-----------------------------------------------------------------------------
static inline
void PMTtoROC( Int_t h_chan,
	       Int_t& crate, Int_t& slot, Int_t& chan )
{
  // Convert location parameters (row, col, chan) of the given PMT
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx
  // In the case of GRINCH/RICH: 
  // crate = GTP; slot = VETROC; chan = PMT. (NINOs are "transparent", in a similar way to the MPDs)

  //Need some form of DB manager !!! or some own DB
  div_t d = div( h_chan , 1); //fManager->GetChanPerSlot() );
  slot = d.quot;
  chan = d.rem;

  d = div( slot, 1); //fManager->GetSlotPerCrate() );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{ //Need some form of DB manager !!! or some own DB
  return chan;// +
  //fManager->GetChanPerSlot()*( slot + fManager->GetSlotPerCrate()*crate );
}

//-----------------------------------------------------------------------------
Int_t SBSSimDecoder::PMTfromROC( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Return index of digitized strip correspomding to hardware channel
  // (crate,slot,chan)

  if( fPMTMap.empty() )
    return -1;

  PMTMap_t::const_iterator found = fPMTMap.find( MakeROCKey(crate,slot,chan) );
  if( found == fPMTMap.end() )
    return -1;

  return found->second;
}

/*
//-----------------------------------------------------------------------------
MCHitInfo SBSSimDecoder::GetMCHitInfo( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Get MC truth info for the given hardware channel
  
  const char* const here = "SBSSimDecoder::GetMCHitInfo";

  // Int_t istrip = StripFromROC( crate, slot, chan );
  // assert( istrip >= 0 );  // else logic error in caller or bad fPMTMap
  
  // assert( buffer );       // Must still have the event buffer
  // const SBSSimEvent* simEvent = reinterpret_cast<const SBSSimEvent*>(buffer);

  // assert( static_cast<vsiz_t>(istrip) < simEvent->fGEMStrips.size() );
  // const SBSSimEvent::DigiGEMStrip& strip = simEvent->fGEMStrips[istrip];
  // assert( strip.fProj >= 0 && strip.fProj < fManager->GetNReadOut() );

  MCHitInfo mc;
  
  // if(strip.fProj==0 && strip.fPlane==4 && strip.fTime1>50.0)
  //   printf("%f \n", strip.fTime1);
  
  if (TESTBIT(strip.fSigType, kInducedStrip) && !TESTBIT(strip.fSigType, kPrimaryStrip) &&
      !TESTBIT(strip.fSigType, kSecondaryStrip) ){
    mc.fMCTrack = 0;
    mc.fMCPos = fManager->GetPosFromSectorStrip(strip.fProj, strip.fPlane, strip.fSector, strip.fChan);
    mc.fMCTime = strip.fTime1;
    //printf("strip num: %d, proj: %d, type: %d, time: %f\n", strip.fChan, strip.fProj, strip.fSigType, strip.fTime1);
    //cout << "strip = " << strip.fChan << ", time = " << mc.fMCTime << ", pos = " <<  mc.fMCPos << endl;
    return mc;
  }
  
  Double_t nOverlapSignal = 0.;
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1
    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const SBSSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
    assert( c.fID == iclust+1 );
    assert( strip.fPlane == c.fPlane && strip.fSector == c.fSector );
    
    Int_t signalID = -1;
    for (unsigned int ii = 0; ii<fSignalInfo.size(); ii++){
      if (c.fType == fSignalInfo.at(ii).tid && c.fPID == fSignalInfo.at(ii).pid)
	signalID = ii;
    }
    // cout << "Plane " << strip.fPlane << ", proj (x: 0, y: 1) " << strip.fProj 
    //  	 << ": pos[proj] = "  << c.fXProj[strip.fProj] << endl;
    if( signalID >= 0 && c.fSource == kPrimarySource ) {
      if( mc.fMCTrack > 0 ) {
        //this means that there two signal hits overlapping
        //for now I keep the fMCTrack to the first one, by average the fMCPos nad fMCTime
        //Weizhi Xiong
        //assert(manager->GetNSigParticle() > 1); //otherwise should not happen
        
        mc.fMCPos += c.fXProj[strip.fProj];
        mc.fMCTime += c.fTime; 
      }else{
        // Strip contains a contribution from a primary particle hit :)
        mc.fMCTrack = fSignalInfo.at(signalID).tid; 
        mc.fMCPos   = c.fXProj[strip.fProj];
        mc.fMCTime  = c.fTime;
      }
      nOverlapSignal++;
    } else {
      ++mc.fContam;
      if( mc.fMCTrack == 0 ) {
	mc.fMCPos += c.fXProj[strip.fProj];
	//printf("strip num: %d, proj: %d, type: %d, time: %f, clust: %d, time: %f \n", strip.fChan, strip.fProj, strip.fSigType, strip.fTime1, iclust, c.fTime);
      }
    }
  }
  assert( strip.fClusters.GetSize() == 0 || mc.fMCTrack > 0 || mc.fContam > 0 );
  
  if( mc.fMCTrack == 0 ) {
    if( mc.fContam > 1 ) {
      // If only background hits, report the mean position of all those hits
      mc.fMCPos /= static_cast<Double_t>(mc.fContam);
    }
    mc.fMCTime = strip.fTime1;
  }else{
    mc.fMCPos /= nOverlapSignal;
    mc.fMCTime /= nOverlapSignal;
  }
  
  
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1
    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const SBSSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
    assert( c.fID == iclust+1 );
    assert( strip.fPlane == c.fPlane && strip.fSector == c.fSector );
    if( c.fType == kPrimaryType && c.fSource == kPrimarySource ) {
      if( mc.fMCTrack > 0 ) {
	Warning( Here(here), "Event %d: Multiple hits of primary particle "
		 "in plane %d\nShould never happen. Call expert.",
		 simEvent->fEvtID, strip.fPlane );
	continue;
      }
      // Strip contains a contribution from a primary particle hit :)
      mc.fMCTrack = 1;    // currently only ever one primary particle per event
      mc.fMCPos   = c.fXProj[strip.fProj];
      mc.fMCTime  = c.fTime;
    } else {
      ++mc.fContam;
      if( mc.fMCTrack == 0 ) {
	mc.fMCPos += c.fXProj[strip.fProj];
      }
    }
  }
  assert( strip.fClusters.GetSize() == 0 || mc.fMCTrack > 0 || mc.fContam > 0 );

  if( mc.fMCTrack == 0 ) {
    if( mc.fContam > 1 ) {
      // If only background hits, report the mean position of all those hits
      mc.fMCPos /= static_cast<Double_t>(mc.fContam);
    }
    mc.fMCTime = strip.fTime1;
  }

  // cout << "SBSSimDecoder: GetMCHitInfo: strip ADC sample: crate " << crate << " slot " << slot << " chan " << chan<< endl;
  // for( int k=0; k<strip.fNsamp; k++ ) {
  //   cout << strip.fADC[k] << " ";
  // }
  // cout << endl;
  
  
  return mc;
}
*/

//-----------------------------------------------------------------------------
static inline Int_t NumberOfSetBits( UInt_t v )
{
  // Count number of bits set in 32-bit integer. From
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

  v = v - ((v >> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
  return (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
Int_t SBSSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t SBSSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'

  static const char* const here = "SBSSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in SBSSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  //const SBSSimEvent* simEvent = reinterpret_cast<const SBSSimEvent*>(buffer);
  
  Int_t ret = HED_OK;

  if (first_decode || fNeedInit) {
    if( (ret = init_cmap()) != HED_OK )
      return ret;
    if( (ret = init_slotdata(fMap)) != HED_OK)
      return ret;
    first_decode = false;
  }
  if( fDoBench ) fBench->Begin("clearEvent");
  Clear();
  for( int i=0; i<fNSlotClear; i++ )
    crateslot[fSlotClear[i]]->clearEvent();
  if( fDoBench ) fBench->Stop("clearEvent");

  // FIXME: needed?
  evscaler = 0;
  event_length = 0;
  
  event_type = 1;
  event_num =  1.;//simEvent->fEvtID;
  recent_event = event_num;

  // Event weight
  fWeight = 1.;//simEvent->fWeight;
  
  //
  if( fDoBench ) fBench->Begin("physics_decode");
  
  Bool_t newclus;
  Int_t crate, slot, chan;
  // Decode the digitized PMT data.  Populate crateslot array.
  /*  for( vector<SBSSimEvent::PMTHit>::size_type i = 0;
       i < simEvent->fPMTHits.size(); ++i ) {
    const SBSSimEvent::PMTHit& h = simEvent->fPMTHits[i];
    
    PMTtoROC(h.fChannel, crate, slot, chan);
    
    // if( c.fPlane < 0 || c.fPlane >= fManager->GetNChamber() ) {
    //   Error( here, "Illegal plane number = %d in cluster. "
    // 	     "Should never happen. Call expert.", c.fPlane );
    //   simEvent->Print("clust");
    //   return HED_FATAL;
    // }
    newclus = true;
    
    new( (*fMCCherHits)[GetNPMThits()] ) SBSSimPMTHit(h);
    
    for(int j = 0; j<GetNPMTclus(); j++){
      
      SBSSimCherCluster* c = GetPMTclus(j);
      
      if(h.fMCtrackID==c->fMCtrackID){
	//Add a hit to the cluster
	newclus = false;
	c->fHitList.push_back(h.fChannel);
	//if( !fHitList ) fHitList = new TList;
	//fHitList->AddLast( theHit );
	c->fSize+=1;
	Int_t listnewsize = c->fSize;
	c->fXcenter = (c->fXcenter*((Double_t)(listnewsize-1))+h.fXPMT)/((Double_t)listnewsize);
	c->fYcenter = (c->fYcenter*((Double_t)(listnewsize-1))+h.fYPMT)/((Double_t)listnewsize);
	
	c->fXcenter_w = c->fXcenter_w*c->fNpe;
	c->fYcenter_w = c->fYcenter_w*c->fNpe;
	c->fNpe+= h.fNpe;
	c->fXcenter_w+= h.fNpe*h.fXPMT;
	c->fYcenter_w+= h.fNpe*h.fYPMT;
	c->fXcenter_w = c->fXcenter_w/c->fNpe;
	c->fYcenter_w = c->fYcenter_w/c->fNpe;
	
	c->fMeanRisingTime = (c->fMeanRisingTime*((Double_t)(listnewsize-1))+h.fTDCtime[0])/((Double_t)listnewsize);
	c->fMeanFallingTime = (c->fMeanFallingTime*((Double_t)(listnewsize-1))+h.fTDCtime[0])/((Double_t)listnewsize);
	c->fRisingTimeRMS = sqrt((pow(c->fRisingTimeRMS, 2)*((Double_t)(listnewsize-1))+ pow(h.fTDCtime[0], 2))/
			      ((Double_t)listnewsize));
	c->fFallingTimeRMS = sqrt((pow(c->fFallingTimeRMS, 2)*((Double_t)(listnewsize-1))+ pow(h.fTDCtime[1], 2))/
			       ((Double_t)listnewsize));
      }
      
    }
    
    if(newclus){
      new( (*fMCCherClus)[GetNPMTclus()] ) SBSSimCherCluster();
      
      
    }
    }*/
  
  /*
  Bool_t newclus = true;
  //fill the 
  for(int i_ = 0; i_<fEvent->fMCClusterHitID.size(); i_++){
    if(tscd.GetParticleType(ih)==fEvent->fMCClusterHitID[i_].first){
      fEvent->fMCClusterHitID[i_].second.push_back(ih);
      newclus = false;
      break;
    }
  }
  if(newclus){
    Int_t trkID = tscd.GetParticleType(ih);
    std::vector<Short_t> ClusterPMTlist;
    ClusterPMTlist.push_back(ih);
    fEvent->fMCClusterHitID.push_back(make_pair(trkID, ClusterPMTlist));
  }
  */
  
  /*
  // Decode the digitized strip data.  Populate crateslot array.
  for( vector<SBSSimEvent::DigiGEMStrip>::size_type i = 0;
       i < simEvent->fGEMStrips.size(); i++) {
    //cout << "i " << i << endl;
    const SBSSimEvent::DigiGEMStrip& s = simEvent->fGEMStrips[i];
    Int_t crate, slot, chan;
    //cout << "striptoroc: " << endl;
    StripToROC( s.fPlane, s.fSector, s.fProj, s.fChan, crate, slot, chan );
    //cout << "crate = " << crate << ", slot = " << slot << ", chan " << chan << endl;
    //cout << "samples: " << endl;
    for( Int_t k = 0; k < s.fNsamp; k++ ) { 
      Int_t raw = s.fADC[k];
      //cout << raw << " ";
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,raw,raw) == SD_ERR )
	return HED_ERR;
    }
    //cout << endl;
    //cout << "stripmap : " << endl;
    // Build map from ROC address to strip index. This is needed to extract
    // the MC truth info later in the tracking detector decoder via GetMCChanInfo.
#ifndef NDEBUG
    pair<PMTMap_t::const_iterator,bool> ins =
#endif
      fPMTMap.insert( make_pair( MakeROCKey(crate,slot,chan), i ) );
    // cout << "ROC key inserted in strip map " << endl;
    // cout << "ins.second ? " << ins.second << endl;
    assert( ins.second );
  }
  
  // Create lists of two types of tracks:
  // 1) Physics tracks, as generated at the target
  // 2) "Back tracks": hits in any GEM plane from the primary particle

  // Physics tracks. We need to copy them here so we can export them as global
  // variables.
  TClonesArray* tracks = simEvent->fMCTracks;
  assert( tracks );
  for( Int_t i = 0; i < tracks->GetLast()+1; i++ ) {
    SBSSimTrack* trk = static_cast<SBSSimTrack*>(tracks->UncheckedAt(i));
   new( (*fMCTracks)[i] ) SBSSimTrack(*trk);
  }
  assert( GetNMCTracks() > 0 );

  // MC hit data ("clusters") and "back tracks"
  Int_t best_primary = -1, best_primary_plane = fManager->GetNChamber(), primary_sector = -1;
  UInt_t primary_hitbits = 0, ufail = 0, vfail = 0;
  for( vector<SBSSimEvent::GEMCluster>::size_type i = 0;
       i < simEvent->fGEMClust.size(); ++i ) {
    const SBSSimEvent::GEMCluster& c = simEvent->fGEMClust[i];

    if( c.fPlane < 0 || c.fPlane >= fManager->GetNChamber() ) {
      Error( here, "Illegal plane number = %d in cluster. "
	     "Should never happen. Call expert.", c.fPlane );
      simEvent->Print("clust");
      return HED_FATAL;
    }

    // Save hits in the GEMs
    new( (*fMCHits)[GetNMCHits()] ) SBSSimGEMHit(c);

    // Extra bookkeeping for primary tracks, used for making back tracks below
    if( c.fType == kPrimaryType && c.fSource == kPrimarySource ) {
      // Record the primary track's points for access via the SimDecoder interface.
      // Record one point per projection so that we can study residuals.
      Int_t itrack = 1;
      primary_sector = c.fSector;
      MCTrackPoint* upt = // kUPlane changed to kXPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kXPlane,
							  c.fMCpos, c.fP );
      upt->fMCTime = c.fTime;
      MCTrackPoint* vpt =// kVPlane changed to kYPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kYPlane,
							  c.fMCpos, c.fP );
      vpt->fMCTime = c.fTime;
      
      // //debug...
      // cout << "SBSSimDecoder.cxx: Print MC points " << endl;
      // cout << "kXplane ? " << kXPlane << endl;
      // upt->Print("");
      // cout << "kVYlane ? " << kYPlane << endl;
      // vpt->Print("");
      
      // Keep bitpattern of planes crossed by this primary
      SETBIT(primary_hitbits,c.fPlane);

      //cout << "Plane number " << c.fPlane << ", primary hitbits " << primary_hitbits << endl; 
      
      // Save index of the primary particle hit closest to plane 0
      if( c.fPlane < best_primary_plane ) {
	best_primary = i;
	best_primary_plane = c.fPlane;
      }
      // Determine digitization hit inefficiency: Check if this MC hit
      // activated GEM strips in both readout planes
      if( c.fSize[0] == 0 ) {
	SETBIT(ufail, c.fPlane);
	CLRBIT(upt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(upt->fStatus, MCTrackPoint::kDigitized);
      }
      if( c.fSize[1] == 0 ) {
	SETBIT(vfail, c.fPlane);
	CLRBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      }
    }
  }
  
  // Sort fMCPoints by type (u,v) and plane number, then calculate plane-to-plane
  // differences. The following assumes that all points are from the same track
  // (ensured above). If that is no longer so one day, fMCPoints will need to
  // be sorted by track number as well, and the algo below needs to be changed.
  fMCPoints->Sort();
  Double_t mass = 0;
  SBSSimTrack* trk = static_cast<SBSSimTrack*>(fMCTracks->UncheckedAt(0));
  assert(trk);
  if( TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(trk->fPID) )
    mass = particle->Mass();
  else
    Warning( "LoadEvent", "No enrty in PDG database for PID = %d", trk->fPID );

  MCTrackPoint* prev_pt = 0;
  for( Int_t i = 0; i < GetNMCPoints(); ++i ) {
    MCTrackPoint* pt = static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(i) );
    assert(pt);
    if( prev_pt && prev_pt->fType == pt->fType ) {
      assert( pt->fMCTrack == prev_pt->fMCTrack );
      if( prev_pt->fPlane+1 == pt->fPlane ) {
	pt->fDeltaE = TMath::Sqrt(prev_pt->fMCP.Mag2() + mass*mass) -
	  TMath::Sqrt(pt->fMCP.Mag2() + mass*mass);
	pt->fDeflect = prev_pt->fMCP.Angle(pt->fMCP);
	pt->fToF = pt->fMCTime - prev_pt->fMCTime;
      }
    }
    prev_pt = pt;
  }

  // Keep statistics in the MC track
  trk->fNHits = 2*NumberOfSetBits(primary_hitbits);
  trk->fHitBits = primary_hitbits;

  // "Back tracks"
  // Record the apparent track from the primary particle
  // of the signal data here, i.e. type == 1 and source == 0.
  // There is only ever one primary particle per event.
  if( best_primary >= 0 ) {

    Int_t nback = GetNBackTracks();
    assert( nback == 0 );

    SBSSimBackTrack* btr = new( (*fBackTracks)[nback] )
      SBSSimBackTrack(simEvent->fGEMClust[best_primary]);

    //cout << "Backtrack primary hitbits " << primary_hitbits << endl;
    
    btr->SetHitBits(primary_hitbits);
    btr->SetUfailBits(ufail);
    btr->SetVfailBits(vfail);

    // Use the back track to emulate calorimeter hits.
    // Assumptions:
    // - Only tracks crossing all fManager->GetNChamber() GEMs (points in all planes)
    //   make a calorimeter hit. This is a crude model for the trigger.
    // - The track propagates without deflection from the last GEM plane
    //   to the front of the emulated calorimeter.
    // - The measured calorimeter position is independent of the incident
    //   track angle.
    if( fManager->DoCalo() && trk->fNHits == 2*fManager->GetNChamber() ) {
      // Retrieve last MC track point
      assert( GetNMCPoints() == 2*fManager->GetNChamber() );
      MCTrackPoint* pt =
	static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(2*fManager->GetNChamber()-1) );
      assert( pt );
      const TVector3& pos = pt->fMCPoint;
      TVector3 dir = pt->fMCP.Unit();
      if( fManager->GetCaloZ() <= pos.Z() ) {
	Error( here, "Calorimeter z = %lf less than z of last GEM plane = %lf. "
	       "Set correct value with SetCaloZ() or turn off calo emulation.",
	       fManager->GetCaloZ(), pos.Z() );
	return HED_FATAL;
      }
      if( TMath::Abs(dir.Z()) < 1e-6 ) {
	Error( here, "Illegal primary track direction (%lf,%lf,%lf). "
	       "Should never happen. Call expert.", dir.X(), dir.Y(), dir.Z() );
	return HED_ERR;
      }
      dir *= 1.0/dir.Z();  // Make dir a transport vector
      TVector3 hitpos = pos + (fManager->GetCaloZ()-pos.Z()) * dir;

      // Smear the position with the given resolution
      // Assumes z-axis normal to calorimeter plane. Otherwise we would have to
      // get the plane's fXax and fYax
      TVector3 res( gRandom->Gaus(0.0, fManager->GetCaloRes()),
		    gRandom->Gaus(0.0, fManager->GetCaloRes()), 0.0 );
      hitpos += res;

      // Encode the raw hit data for the dummy GEM planes.
      // The actual coordinate transformation to u or v takes place in each
      // plane's Decode() where all the required geometry information is at hand.
      // This bypasses any type of digitization. That should be fine for dummy
      // planes where we want to inject known lab hit positions into the tracking.
      //
      // Because of the way the detector map is layed out at the moment,
      // we place the calorimeter in fake sector 31, so the data are in two slots
      // (for u and v coordinates, respectively) in the ROC immediately
      // following the GEM trackers for sector 30. In each slot, channels
      // 0-29 correspond to the sector of the MC track sector (should always be 0
      // if mapping sectors. Each "hit" corresponds to one measured position.
      // Currently, there is only ever one hit per channel since there is only
      // one MC track. The hit's raw data are hitpos.X(), the data, hitpos.Y(),
      // each a Float_t value interpreted as Int_t.
      assert( primary_sector == 0 );

      union FloatIntUnion {
	Float_t f;
	Int_t   i;
      } datx, daty;
      datx.f = static_cast<Float_t>(hitpos.X());
      daty.f = static_cast<Float_t>(hitpos.Y());

      Int_t crate, slot, chan;
      //StripToROC( 0, fManager->GetNSector(), kUPlane, primary_sector, crate, slot, chan );
      StripToROC( 0, fManager->GetNSector(), kXPlane, primary_sector, crate, slot, chan );
       if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	return HED_ERR;
      //StripToROC( 0, fManager->GetNSector(), kVPlane, primary_sector, crate, slot, chan );
      StripToROC( 0, fManager->GetNSector(), kYPlane, primary_sector, crate, slot, chan );
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	return HED_ERR;
    }
  }

  // DEBUG:
  //cout << "SimDecoder: nTracks = " << GetNMCTracks() << endl;
  //fMCTracks.Print();
  */

  return HED_OK;
}
/*
//-----------------------------------------------------------------------------
SBSSimPMTHit::SBSSimPMTHit( const SBSSimEvent::PMTHit& h )
  : fID(h.fID), 
    fSource(h.fSource), fType(h.fType), 
    fMCtrackPID(h.fMCtrackPID), fOrigVolFlag(h.fOrigVolFlag), 
    fXPMT(h.fXPMT), fYPMT(h.fYPMT), 
    fNpe(h.fNpe), fTime(h.fTime), 
    fTDCtime1(h.fTDCtime[0]), fTDCtime2(h.fTDCtime[1]),
    fDetID(h.fDetID), fChannel(h.fChannel), 
    fPMTrow(h.fPMTrow), fPMTcol(h.fPMTcol), 
    fTDC1(h.fTDC[0]), fTDC2(h.fTDC[1])
{
  // Construct hit from cluster
}
*/
//-----------------------------------------------------------------------------
void SBSSimPMTHit::Print( const Option_t* ) const
{
  // Print SBSSimGEMHit info
}

ClassImp(SBSSimPMTHit)
ClassImp(SBSSimCherCluster)
