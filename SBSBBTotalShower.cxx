///////////////////////////////////////////////////////////////////////////////
//                                                                           
// THaBBTotalShower                                                            
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

#include "THaBBTotalShower.h"
#include "THaBBShower.h"

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

ClassImp(THaBBTotalShower)

//_____________________________________________________________________________
THaBBTotalShower::THaBBTotalShower( const char* name, const char* description,
                                   THaApparatus* apparatus ) :
THaPidDetector(name,description,apparatus), 
fShower(NULL), fPreShower(NULL), fMaxDx(0.0), fMaxDy(0.0), fE(NULL)
{
    // Constructor. With this method, the subdetectors are created using
    // this detector's prefix followed by "sh" and "ps", respectively,
    // and variable names like "L.ts.sh.nhits".

    DEBUG_HALL_A_ANALYZER_DEBUGER_INIT;
    DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER;

    Setup( name, "sh", "ps", description, apparatus, true );
}


//_____________________________________________________________________________
THaBBTotalShower::THaBBTotalShower( const char* name, 
                                   const char* shower_name,
                                   const char* preshower_name,
                                   const char* description,
                                   THaApparatus* apparatus ) :
THaPidDetector(name,description,apparatus),
fShower(NULL), fPreShower(NULL), fMaxDx(0.0), fMaxDy(0.0), fE(NULL)
{
    // Constructor. With this method, the subdetectors are created using
    // the given names 'shower_name' and 'preshower_name', and variable 
    // names like "L.sh.nhits".

    Setup( name, shower_name, preshower_name, description, apparatus, false );
}

//_____________________________________________________________________________
void THaBBTotalShower::Setup( const char* name,
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

    fShower = new THaBBShower( sname, desc, apparatus );
    if( !fShower || fShower->IsZombie() ) {
        MakeZombie();
        goto exit;
    }

    if( subnames )
        strcpy( subname+nlen+1, preshower_name );
    else
        sname = preshower_name;
    strcpy( desc+dlen, " preshower subdetector" );

    fPreShower = new THaBBShower( sname, desc, apparatus );
    if( !fPreShower && fPreShower->IsZombie() ) {
        MakeZombie();
        goto exit;
    }

exit:
    delete [] subname;
    delete [] desc;
    return;
}

//_____________________________________________________________________________
void THaBBTotalShower::ClearEvent() {

    for (Int_t i=0;i<kMaxNClust;i++) {
        fE[i] = 0.0;
        fX[i] = kBig;
        fY[i] = kBig;
        fID[i] = 0;
    }
}

//_____________________________________________________________________________
THaBBTotalShower::~THaBBTotalShower()
{
    // Destructor. Remove variables from global list.
    if( fIsSetup )
        RemoveVariables();

    delete [] fE; fE = 0;
    delete [] fX; fX = 0;
    delete [] fY; fY = 0;
    delete [] fID; fID = 0;

    delete fPreShower;
    delete fShower;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus THaBBTotalShower::Init( const TDatime& run_time )
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
Int_t THaBBTotalShower::ReadDatabase( const TDatime& date )
{
    // Read this detector's parameters from the database file 'fi'.
    // This function is called by THaDetectorBase::Init() once at the
    // beginning of the analysis.
    // 'date' contains the date/time of the run being analyzed.

    const int LEN = 100;
    char line[LEN];

    FILE* fi = OpenFile( date );
    if( !fi ) return kFileError;

    fgets ( line, LEN, fi ); fgets ( line, LEN, fi );          
    fscanf ( fi, "%f%f", &fMaxDx, &fMaxDy );  // Max diff of shower centers

    fIsInit = true;
    fclose(fi);

    RemoveVariables();
    if (fE) {
        delete [] fE;   fE = 0;
        delete [] fX;   fX = 0;
        delete [] fY;   fY = 0;
        delete [] fID;  fID = 0;
    }

    fE = new Float_t[kMaxNClust];
    fX = new Float_t[kMaxNClust];
    fY = new Float_t[kMaxNClust];
    fID = new Int_t[kMaxNClust];

    return kOK;
}

//_____________________________________________________________________________
Int_t THaBBTotalShower::DefineVariables( EMode mode )
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
        { "id", "ID of Psh&Sh coincidence (1==good)", "fID" },
        { 0 }
    };
    return DefineVarsFromList( vars, mode );
}

//_____________________________________________________________________________
Int_t THaBBTotalShower::Decode( const THaEvData& evdata )
{
    // Decode total shower detector. Calls Decode() of fPreShower and fShower.
    // Return the return value of fShower->Decode().

    if( !IsOK() ) 
        return -1;

    ClearEvent();

    fPreShower->Decode( evdata );
    return fShower->Decode( evdata );
}

//_____________________________________________________________________________
Int_t THaBBTotalShower::CoarseProcess(TClonesArray& tracks )
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

    fPreShower->CoarseProcess(tracks );
    fShower->CoarseProcess( tracks );

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
}
//_____________________________________________________________________________
Int_t THaBBTotalShower::FineProcess( TClonesArray& tracks )
{
    // Fine processing. 
    // Call fPreShower->FineProcess() and fShower->FineProcess() in turn.
    // Return return value of fShower->FineProcess().

    if( !IsOK() )
        return -1;

    fPreShower->FineProcess( tracks );
    return fShower->FineProcess( tracks );
}

//_____________________________________________________________________________
void THaBBTotalShower::SetApparatus( THaApparatus* app )
{
    // Set the apparatus of this detector as well as the subdetectors

    THaPidDetector::SetApparatus( app );
    fShower->SetApparatus( app );
    fPreShower->SetApparatus( app );
    return;
}

//_____________________________________________________________________________

void THaBBTotalShower::LoadMCHitAt( Double_t x, Double_t y, Double_t E )
{
    ClearEvent();
    fNclust = 0;
    fE[fNclust] = E;
    fX[fNclust] = x;
    fY[fNclust] = y;
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
