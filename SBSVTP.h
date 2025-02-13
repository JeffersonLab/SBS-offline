#ifndef ROOT_SBSVTP
#define ROOT_SBSVTP

#include "THaDetector.h"

// VTP Data handler

struct VTPCluster {
  UInt_t fDet;           // Detector ID
  vector<Int_t> fX;      // cluster x coord, only for seed tower
  vector<Int_t> fY;      // cluster y coord
  vector<Int_t> fE;      // cluster energy
  vector<Int_t> fTime;   // cluster time
  vector<Int_t> fSize;   // cluster size (numbe of hits)
  void clear()
  {
    fDet = 0;
    fX.clear(); fY.clear(); fE.clear(); fTime.clear(); fSize.clear();
  }
};

class SBSVTP : public THaDetector {
  public:
    SBSVTP(const char* name, const char* description = "", THaApparatus* apparatus = nullptr);
    virtual ~SBSVTP();
    
    virtual Int_t Decode( const THaEvData& );
    virtual void  Clear( Option_t* opt="" );
    virtual Int_t DefineVariables( EMode mode );

    const VTPCluster& GetClusters() { return fVTPClusters; }    

  protected:
    VTPCluster fVTPClusters;

  ClassDef(SBSVTP,0);
};

#endif

