///////////////////////////////////////////////////////////////////////////////
//                                                                           
// SBSBBTotalShower                                                            
//                                                                           
// A total shower counter, consisting of a shower and a preshower.           
// Calculates the total energy deposited in Shower+Preshower.           
//                                                                           
// This class is transplanted from AGen offline analysis                     
// http://hallaweb.jlab.org/experiment/E02-013/offlana.html                  
//                                                                           
///////////////////////////////////////////////////////////////////////////////
//
// Notice:  Scintillator planes working with shower detector uses class
// THaScintillator in standard analyzer. They could be actived in scripts
// by include following lines in your replay scripts:
//
//  BB->AddDetector( new THaScintillator( "s", "scintillator" ) );
//  BB->AddDetector( new THaScintillator("sum","BigBite Total Sum") );
//  BB->AddDetector( new THaScintillator("psum","BigBite PS Sum") ); 
//                                                                           
///////////////////////////////////////////////////////////////////////////////

#include "THaApparatus.h"
#include "SBSBBTotalShower.h"
#include "SBSBBShower.h"
#include "SBSCalorimeter.h"

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
THaShower(name,description,apparatus), 
  fShower(NULL), fPreShower(NULL)//, //fMaxDx(0.0), fMaxDy(0.0), 
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
  THaShower(name,description,apparatus),
  fShower(NULL), fPreShower(NULL)//, //fMaxDx(0.0), fMaxDy(0.0), fE(NULL)
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
void SBSBBTotalShower::ClearEvent() {
  fShower->Clear();
  fPreShower->Clear();
  
  fNclust = 0;
  fE = 0.0;
  fX = 0.0;
  fY = 0.0;
  /*
    for (Int_t i=0;i<kMaxNClust;i++) {
        fE[i] = 0.0;
        fX[i] = kBig;
        fY[i] = kBig;
        fID[i] = 0;
    }
  */
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
    if( (status = THaPidDetector::Init( run_time )) ||
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

    //I don't know if this is right, but the other way was wrong
    //The variables defined here blow everything up

    RVarDef vars[] = {
        { "e",  "Energy (MeV) of largest cluster",    "fE" },
        { "x",  "Energy (MeV) of largest cluster",    "fX" },
        { "y",  "Energy (MeV) of largest cluster",    "fY" },
        //{ "id", "ID of Psh&Sh coincidence (1==good)", "fID" },
        { 0 }
    };
    return DefineVarsFromList( vars, mode );
}

//_____________________________________________________________________________
Int_t SBSBBTotalShower::Decode( const THaEvData& evdata )
{
    // Decode total shower detector. Calls Decode() of fPreShower and fShower.
    // Return the return value of fShower->Decode().

    if( !IsOK() ) 
        return -1;

    ClearEvent();
    
    fShower->Decode( evdata );
    return fPreShower->Decode( evdata );
}

//_____________________________________________________________________________
Int_t SBSBBTotalShower::CoarseProcess(TClonesArray& tracks )
{
  // Reconstruct Clusters in shower and preshower detectors.
  // Then compute total shower energy and cluster ID.
  //
  // Valid fIDs:
  //       1   Matching clusters found.
  //       0   Clusters found, but separated too far
  //      -1   No cluster found in either shower or preshower or both
  //
  if( !IsOK() ) 
    return -1;

  // EPAF (2020/05/11): Different strategy
  // Process shower first (since it is *much* cleaner)
  fShower->CoarseProcess( tracks );
  
  if(fShower->GetNclust()){
    // Then set constraints on preshower 
    int rowmin = fPSSHmatchmapX.at(fShower->GetRowMax()).first;
    int rowmax = fPSSHmatchmapX.at(fShower->GetRowMax()).second;
    int colmin = fPSSHmatchmapY.at(fShower->GetColMax()).first;
    int colmax = fPSSHmatchmapY.at(fShower->GetColMax()).second;
    //fShower->GetY();
    
    fPreShower->SetSearchRegion(rowmin, rowmax, colmin, colmax);
    // Finally process shower
    fPreShower->CoarseProcess(tracks );
    
    if(fShower->GetNclust()){
      fNclust = 1;
      /*
      cout << " fShower->GetE() " << fShower->GetE() 
	   << " fPreShower->GetE() " << fPreShower->GetE() << endl;
	   
      fE = fShower->GetE()+fPreShower->GetE();
      
      double w2sh = pow(fShower->GetE()*fShower->GetXSize()/fShower->GetNRows(), 2);
      double w2ps = pow(fPreShower->GetE()*fPreShower->GetXSize()/fPreShower->GetNRows(), 2);

      fX = (fShower->GetX()*w2sh+fPreShower->GetX()*w2ps)/(w2sh+w2ps);
      
      w2sh = pow(fShower->GetE()*fShower->GetYSize()/fShower->GetNCols(), 2);
      w2ps = pow(fPreShower->GetE()*fPreShower->GetYSize()/fPreShower->GetNCols(), 2);

      fY = (fShower->GetY()*w2sh+fPreShower->GetY()*w2ps)/(w2sh+w2ps);
      
      cout << " fNclust " << fNclust << " fE " << fE << " fX " << fX << " fY " << fY << endl;
      */
      return 1;
    }else{
      return -1;
    }
  }else{
    return -1;
  }
  
  /*
    fNclust = 0;
    
    //  if( ( fShower->GetNclust()+fPreShower->GetNclust() )> kMaxNClust ) cerr << "Error:  Too many clusters ( " << fShower->GetNclust() << " shower and " << fPreShower->GetNclust() << " preshower ).  We are restricted to " << kMaxNClust << endl;
    for (int sh=0;sh<fShower->GetNclust();sh++) {
        Bool_t keep = false;
        Int_t ps = 0;
        while ((ps < fPreShower->GetNclust()) && (fNclust < kMaxNClust)) {
            double dx = 
                fShower->GetClust(sh)->GetX()-fPreShower->GetClust(ps)->GetX();
            double dy =
                fShower->GetClust(sh)->GetY()-fPreShower->GetClust(ps)->GetY();
            if (TMath::Abs(dx)<fMaxDx && TMath::Abs(dy)<fMaxDy && (fNclust < kMaxNClust)) {
                keep = true;
                fE[fNclust] = fShower->GetClust(sh)->GetE() + fShower->GetClust(ps)->GetE();
                fX[fNclust] = fShower->GetClust(sh)->GetX();
                fY[fNclust] = fShower->GetClust(sh)->GetY();
                fID[fNclust] = 1;
                fNclust++;
            }
            ps++;
        }
        // This throws away if they cluster positions are not within some range.  Let's not do that for now
        //    if (!keep) 
        //fShower->RemoveCluster(sh--);
    }

    for (int ps=0;ps<fPreShower->GetNclust();ps++) {
        Bool_t keep = false;
        Int_t sh = 0;
        while (!keep && sh < fShower->GetNclust()) {
            double dx = 
                fShower->GetClust(sh)->GetX()-fPreShower->GetClust(ps)->GetX();
            double dy =
                fShower->GetClust(sh)->GetY()-fPreShower->GetClust(ps)->GetY();
            if (TMath::Abs(dx)<fMaxDx && TMath::Abs(dy)<fMaxDy) 
                keep = true;
            sh++;
        }
        // This throws away if they cluster positions are not within some range.  Let's not do that for now
        //if (!keep) 
        //      fPreShower->RemoveCluster(ps--);
    } 


    return 0;
  */
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
    
    if(fNclust){
      //cout << " fShower->GetE() " << fShower->GetE() << " fPreShower->GetE() " << fPreShower->GetE() << endl;
      
      fE = fShower->GetE()+fPreShower->GetE();
      
      double w2sh = pow(fShower->GetE()*fShower->GetXSize()/fShower->GetNRows(), 2);
      double w2ps = pow(fPreShower->GetE()*fPreShower->GetXSize()/fPreShower->GetNRows(), 2);
      
      fX = (fShower->GetX()*w2sh+fPreShower->GetX()*w2ps)/(w2sh+w2ps);
      
      w2sh = pow(fShower->GetE()*fShower->GetYSize()/fShower->GetNCols(), 2);
      w2ps = pow(fPreShower->GetE()*fPreShower->GetYSize()/fPreShower->GetNCols(), 2);
      
      fY = (fShower->GetY()*w2sh+fPreShower->GetY()*w2ps)/(w2sh+w2ps);
      
      //cout << " fNclust " << fNclust << " fE " << fE << " fX " << fX << " fY " << fY << endl;  
    }
    return 0;
}

//_____________________________________________________________________________
void SBSBBTotalShower::SetApparatus( THaApparatus* app )
{
    // Set the apparatus of this detector as well as the subdetectors

    THaPidDetector::SetApparatus( app );
    fShower->SetApparatus( app );
    fPreShower->SetApparatus( app );
    return;
}

//_____________________________________________________________________________

void SBSBBTotalShower::LoadMCHitAt( Double_t x, Double_t y, Double_t E )
{
    ClearEvent();
    fNclust = 0;
    fE = E;
    fX = x;
    fY = y;
    //fE[fNclust] = E;
    //fX[fNclust] = x;
    //fY[fNclust] = y;
    fNclust++;
    if( fShower )
    {
        fShower->LoadMCHitAt( x, y, E );
    }
    if( fPreShower )
    {
        fPreShower->LoadMCHitAt( x, y, E );
    }

}

///////////////////////////////////////////////////////////////////////////////
