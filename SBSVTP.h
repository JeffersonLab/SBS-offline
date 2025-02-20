#ifndef ROOT_SBSVTP
#define ROOT_SBSVTP

#include "THaNonTrackingDetector.h"
#include <vector>

// VTP Data handler

struct VTPCluster {
  UInt_t fDet;           // Detector ID
  std::vector<UInt_t> fX;      // cluster x coord, only for seed tower
  std::vector<UInt_t> fY;      // cluster y coord
  std::vector<UInt_t> fE;      // cluster energy
  std::vector<UInt_t> fTime;   // cluster time
  std::vector<UInt_t> fSize;   // cluster size (numbe of hits)
  void clear()
  {
    fDet = 0;
    fX.clear(); fY.clear(); fE.clear(); fTime.clear(); fSize.clear();
  }
};

class SBSVTP : public THaNonTrackingDetector {
public:
  SBSVTP(const char* name, const char* description = "", THaApparatus* apparatus = nullptr);
  virtual ~SBSVTP();
    
  virtual Int_t Decode( const THaEvData& );
  virtual Int_t CoarseProcess( TClonesArray& tracks);
  virtual Int_t FineProcess( TClonesArray& tracks);
  virtual void  Clear( Option_t* opt="" );
  virtual Int_t ReadDatabase( const TDatime& date );
  virtual Int_t DefineVariables( EMode mode );

  const VTPCluster& GetClusters() { return fVTPClusters; }    

protected:
  Int_t fVTPErrorFlag;
  VTPCluster fVTPClusters;

  ClassDef(SBSVTP,0);
};

#endif

