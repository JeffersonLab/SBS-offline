///////////////////////////////////////////////////////////////////////////////
//                                                                           
// SBSBBTotalShower                                                            
//                                                                           
// A total shower counter, consisting of a shower and a preshower.           
// Calculates the total energy deposited in Shower+Preshower.           
//                                                                           
//                                                                           

#include "THaApparatus.h"
#include "THaPidDetector.h"
#include "SBSBBTotalShower.h"
#include "SBSBBShower.h"

#include "VarType.h"
#include "VarDef.h"
#include "TMath.h"
#include <cstring>
#include <iostream>
#include <iomanip>

using namespace std;

class THaEvData;
class TClonesArray;
class TDatime;

ClassImp(SBSBBTotalShower)

//_____________________________________________________________________________
SBSBBTotalShower::SBSBBTotalShower( const char* name, const char* description,
                                   THaApparatus* apparatus ) :
SBSCalorimeter(name,description,apparatus), 
  fShower(NULL), fPreShower(NULL), fMaxDx(0.4), fMaxDy(0.4)
//fE(0.0), fX(0.0), fY(0.0)//, fID(NULL)
{
    // Constructor. With this method, the subdetectors are created using
    // this detector's prefix followed by "sh" and "ps", respectively,
    // and variable names like "L.ts.sh.nhits".

    DEBUG_HALL_A_ANALYZER_DEBUGER_INIT;
    DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER;

    Setup( name, "sh", "ps", description, apparatus, true );
}


//_____________________________________________________________________________
SBSBBTotalShower::SBSBBTotalShower( const char* name, 
				    const char* shower_name,
                                   const char* preshower_name,
                                   const char* description,
                                   THaApparatus* apparatus ) :
  SBSCalorimeter(name,description,apparatus),
  fShower(NULL), fPreShower(NULL),fMaxDx(0.4), fMaxDy(0.4)
  //fE(0.0), fX(0.0), fY(0.0)
{
    // Constructor. With this method, the subdetectors are created using
    // the given names 'shower_name' and 'preshower_name', and variable 
    // names like "L.sh.nhits".

    Setup( name, shower_name, preshower_name, description, apparatus, false );
}

//_____________________________________________________________________________
void SBSBBTotalShower::Setup( const char* name,
                             const char* shower_name,
                             const char* preshower_name,
                             const char* description,
                             THaApparatus* apparatus,
                             bool subnames )
{
    DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER;
    DEBUG_HALL_A_ANALYZER_DEBUGER_INIT;

    // Set up the total shower counter. Called by constructor.

    static const char* const here = "Setup()";
    static const char* const message = 
        "Must construct %s detector with valid name! Object construction failed.";

    // Base class constructor failed?
    if( IsZombie()) return;

    size_t sh, ps;
    if( !shower_name || (sh = strlen(shower_name)) == 0 ) {
        Error( Here(here), message, "shower" );
        MakeZombie();
        return;
    }
    if( !preshower_name || (ps = strlen(preshower_name)) == 0 ) {
        Error( Here(here), message, "preshower" );
        MakeZombie();
        return;
    }

    size_t nlen = strlen(name);
    size_t slen = TMath::Max(ps,sh);
    size_t len = slen;
    if( subnames )
        len += nlen+1;
    char* subname = new char[ len+1 ];
    const char* sname;
    if( subnames ) {
        strcpy( subname, name );
        strcat( subname, "." );
        strcat( subname, shower_name );
        sname = subname;
    } else 
        sname = shower_name;

    char* desc = new char[ 50+strlen(description) ];
    if( description && *description )
        strcpy( desc, description );
    else {
        strcpy( desc, "Total shower counter" );
        SetTitle( desc );
    }
    size_t dlen = strlen(desc);
    strcat( desc, " shower subdetector" );

    fShower = new SBSBBShower( sname, desc, apparatus );
    //fShower = new SBSCalorimeter( sname, desc, apparatus );
    if( !fShower || fShower->IsZombie() ) {
        MakeZombie();
        goto exit;
    }
    //GetApparatus()->AddDetector(fShower);
    //if(GetApparatus())cout << GetApparatus()->GetName() << endl;
    //if(fShower->GetApparatus())cout << fShower->GetApparatus()->GetName() << endl;

    if( subnames )
        strcpy( subname+nlen+1, preshower_name );
    else
        sname = preshower_name;
    strcpy( desc+dlen, " preshower subdetector" );

    fPreShower = new SBSBBShower( sname, desc, apparatus );
    //fPreShower = new SBSCalorimeter( sname, desc, apparatus );
    if( !fPreShower && fPreShower->IsZombie() ) {
        MakeZombie();
        goto exit;
    }
    //GetApparatus()->AddDetector(fShower);
    
exit:
    delete [] subname;
    delete [] desc;
    return;
}
//_____________________________________________________________________________
Int_t SBSBBTotalShower::Decode( const THaEvData& evdata )
{
  fShower->Decode(evdata);
  fPreShower->Decode(evdata);
  return 0;
}
//_____________________________________________________________________________
void SBSBBTotalShower::ClearEvent() {
  fShower->ClearEvent();
  fPreShower->ClearEvent();
  
}

//_____________________________________________________________________________
SBSBBTotalShower::~SBSBBTotalShower()
{
    // Destructor. Remove variables from global list.
    if( fIsSetup )
        RemoveVariables();
    
    /*
    delete [] fE; fE = 0;
    delete [] fX; fX = 0;
    delete [] fY; fY = 0;
    delete [] fID; fID = 0;
    */
    fPSSHmatchmapX.clear();
    fPSSHmatchmapY.clear();
    
    //delete fPSSHmatchmapX;
    //delete fPSSHmatchmapY;
    
    delete fPreShower;
    delete fShower;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus SBSBBTotalShower::Init( const TDatime& run_time )
{
    // Set up total shower counter. This function 
    // reads the total shower parameters from a local database file 
    // (e.g. "db_ts.dat"), then calls Init() for
    // the shower and preshower subdetectors, respectively.

    if( IsZombie() || !fShower || !fPreShower )
        return fStatus = kInitError;

    EStatus status;
    if( (status = SBSCalorimeter::Init( run_time )) ||
        (status = fShower->Init( run_time )) ||
        (status = fPreShower->Init( run_time )) )
        return fStatus = status;

    return fStatus;
}

//_____________________________________________________________________________
Int_t SBSBBTotalShower::ReadDatabase( const TDatime& date )
{
  // Read this detector's parameters from the database file 'fi'.
  // This function is called by THaDetectorBase::Init() once at the
  // beginning of the analysis.
  // 'date' contains the date/time of the run being analyzed.
  
  static const char* const here = "ReadDatabase()";
  
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  
  /*
  // Read fOrigin and fSize (required!) // ou pas
  Int_t err = ReadGeometry( file, date, true );
  if( err ) {
    //cout << " readgeo err" << endl;
    fclose(file);
    return err;
  }
  */
  
  string components_names;
  std::vector<Int_t> pssh_matchmap_x;
  std::vector<Int_t> pssh_matchmap_y;
  
  DBRequest config_request[] = {
    { "components",   &components_names,   kString, 0, 1},
    { "pssh_matchmap_x",   &pssh_matchmap_x,   kIntV, 0, 1},
    { "pssh_matchmap_y",   &pssh_matchmap_y,   kIntV, 0, 1},
    { 0 }
  };
  
  Int_t err = LoadDB( file, date, config_request, fPrefix );

  if( err ) {
    return kInitError;
  }
  
  for(int i = 0; i<pssh_matchmap_x.size(); i+=3){
    fPSSHmatchmapX[ pssh_matchmap_x[i] ] = 
      std::make_pair(pssh_matchmap_x[i+1], pssh_matchmap_x[i+2]);
  }

  for(int i = 0; i<pssh_matchmap_y.size(); i+=3){
    fPSSHmatchmapY[ pssh_matchmap_y[i] ] = 
      std::make_pair(pssh_matchmap_y[i+1], pssh_matchmap_y[i+2]);
  }
  
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSBBTotalShower::DefineVariables( EMode mode )
{
    // Initialize global variables and lookup table for decoder

    if( mode == kDefine && fIsSetup ) return kOK;
    fIsSetup = ( mode == kDefine );

    // Register global variables

    return 0;
}

//_____________________________________________________________________________
Int_t SBSBBTotalShower::CoarseProcess(TClonesArray& tracks )
{
  //   cout << "SBSBBTotalShower::CoarseProcess " << endl;

  fShower->CoarseProcess(tracks);
  fPreShower->CoarseProcess(tracks);
  //

  fShower->MakeGoodBlocks();
  fPreShower->MakeGoodBlocks();

  //
  fShower->FindClusters();
  // match blocks hit in Preshower to clusters in the  Shower
  std::vector<SBSCalorimeterCluster*> ShowerClusters = fShower->GetClusters();
  std::vector<SBSBlockSet> PreShowerBlockSet = fPreShower->GetBlockSet();
  Int_t PreShower_Nclus= 0;
  for (Int_t nc=0;nc<ShowerClusters.size();nc++) {
    Float_t xsh = ShowerClusters[nc]->GetX();
    Float_t ysh = ShowerClusters[nc]->GetY();
    Bool_t AddToPreShowerCluster = kFALSE;
    for (Int_t nps=0;nps<PreShowerBlockSet.size();nps++) {
      if (!PreShowerBlockSet[nps].InCluster) {
	SBSElement* psblk = fPreShower->GetElement(PreShowerBlockSet[nps].id);
	     Float_t xps =  PreShowerBlockSet[nps].x;
	     Float_t yps =  PreShowerBlockSet[nps].y;
      Bool_t MatchCriterion = abs(xsh-xps) < fMaxDx;
      if (MatchCriterion) {
	PreShowerBlockSet[nps].InCluster = kTRUE;
	if (!AddToPreShowerCluster) {
	  fPreShower->MakeCluster(PreShowerBlockSet.size(),psblk);
	  AddToPreShowerCluster = kTRUE;
	  PreShower_Nclus++;
	} else {
	  fPreShower->AddToCluster(PreShower_Nclus-1,psblk);
	}
      }
      }
    }
    if (!AddToPreShowerCluster && PreShowerBlockSet.size()>0) fPreShower->MakeCluster(PreShowerBlockSet.size()); // If preshower not matched to shower, make preshower cluster with mult = 0
  }
    //
    fPreShower->MakeMainCluster();
  //
  return 0;
}
//_____________________________________________________________________________
Int_t SBSBBTotalShower::FineProcess( TClonesArray& tracks )
{
    // Fine processing. 
    // Call fPreShower->FineProcess() and fShower->FineProcess() in turn.
    // Return return value of fShower->FineProcess().

    if( !IsOK() )
        return -1;

    fPreShower->FineProcess( tracks );
    fShower->FineProcess( tracks );
    
    return 0;
}

//_____________________________________________________________________________
void SBSBBTotalShower::SetApparatus( THaApparatus* app )
{
    // Set the apparatus of this detector as well as the subdetectors

    SBSCalorimeter::SetApparatus( app );
    fShower->SetApparatus( app );
    fPreShower->SetApparatus( app );
    return;
}

//_____________________________________________________________________________

void SBSBBTotalShower::LoadMCHitAt( Double_t x, Double_t y, Double_t E )
{
    ClearEvent();
    /*
  fNclust = 0;
    fE = float(E);
    fX = x;
    fY = y;
    fNclust++;
    if( fShower )
    {
        fShower->LoadMCHitAt( x, y, E );
    }
    if( fPreShower )
    {
        fPreShower->LoadMCHitAt( x, y, E );
    }

    */
}

///////////////////////////////////////////////////////////////////////////////
