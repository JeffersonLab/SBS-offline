#ifndef SBSGEMSTAND_H
#define SBSGEMSTAND_H
#include <vector>
#include <THaTrackingDetector.h>

class THaRunBase;
class THaApparatus;
class THaEvData;
class SBSGEMPlane;
class THaCrateMap;

class SBSGEMStand : public THaTrackingDetector {
    public:
        SBSGEMStand( const char *name, const char *description = "",
                THaApparatus *app = 0 );

        virtual ~SBSGEMStand();

        virtual void    Clear( Option_t* opt="" );
        virtual Int_t   Decode( const THaEvData& );
        virtual EStatus Init( const TDatime& date );

        virtual Int_t   ReadDatabase( const TDatime& date );

        virtual Int_t   CoarseTrack( TClonesArray& tracks );
        virtual Int_t   FineTrack( TClonesArray& tracks );
        virtual Int_t   DefineVariables( EMode mode = kDefine );
        virtual void    Print(const Option_t* opt) const;
        virtual void    SetDebug( Int_t level );

        virtual Int_t   Begin( THaRunBase* r=0 );
        virtual Int_t   End( THaRunBase* r=0 );

    private:
        std::vector <SBSGEMPlane *> fPlanes;
	
	bool fIsMC;
	
        THaCrateMap *fCrateMap;
        ClassDef(SBSGEMStand ,0);
};

#endif
