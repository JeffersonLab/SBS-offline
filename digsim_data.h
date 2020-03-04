#ifndef digsim_data_h
#define digsim_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
class TTree;


namespace SBSDigSim {
  //dummy class to get rid of "unused class rule"
  class digsim_data {
  public :
    digsim_data(){};
    virtual ~digsim_data(){};
  };

  // Purely virtual data structure for detector data in a G4SBS ROOT tree.
  struct VDetData_t {
    VDetData_t() {};
    virtual ~VDetData_t(){};
    // All sub-classes *must* implement a concrete instance of this
    virtual bool SetupBranches(TTree *t, const char* prefix) = 0;
    template<typename T>
    int SetupBranch(TTree* tree,const char* prefix,const char* varname,T &var);
  };
  
  //To read the "MC truth" info for hits from PMT detectors
  struct PMTSimHit_t : public VDetData_t {
    UInt_t                nhits;   // number of hits
    std::vector<short>   *src;     // source (sig, bkgd)
    std::vector<short>   *trid;    // track ID
    std::vector<int>     *pid;     // Particle ID
    std::vector<short>   *chan;    // channel (PMT)
    std::vector<double>  *edep;    // energy deposit in element (for non Ckov dets)
    std::vector<int>     *npe;     // number of photoelectrons
    std::vector<double>  *time;    // hit time
    std::vector<double>  *t_lead;  // lead time (for dets w/ TDCs)
    std::vector<double>  *t_trail; // trail time (for dets w/ TDCs)
    bool fReadEdep;                // flag to read edep. Set to true by default
    bool fReadTimes;               // flag to read times. Set to true by default
    // these flags are set in constructor - just below
    PMTSimHit_t(bool readTimes = true, bool readEdep = true) : nhits(0),
      src(0), trid(0), pid(0), 
      chan(0), edep(0), npe(0), 
      time(0), t_lead(0), t_trail(0), 
      fReadEdep(readEdep), fReadTimes(readTimes)
      {}
    virtual ~PMTSimHit_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  //To read the "MC truth" info for hits from PMT detectors
  struct GEMSimHit_t : public VDetData_t {
    UInt_t                nhits;   // number of hits
    std::vector<short>   *src;     // source (sig, bkgd)
    std::vector<short>   *trid;    // track ID
    std::vector<int>     *pid;     // Particle ID
    std::vector<short>   *plane;   // GEM plane
    std::vector<short>   *module;  // GEM module
    std::vector<double>  *edep;    // energy deposit
    std::vector<double>  *time;    // hit time
    std::vector<double>  *xpos;    // true MC x position
    std::vector<double>  *ypos;    // true MC y position
    // momentum of particle: x, y, z
    std::vector<double>  *px;      
    std::vector<double>  *py;
    std::vector<double>  *pz;
    // hit extension (in number of strips) in both projections
    std::vector<short>   *sizex;
    std::vector<short>   *sizey;
    // hit extension (in number of strips) in both projections
    std::vector<short>   *startx;
    std::vector<short>   *starty;
    GEMSimHit_t() : nhits(0),
      src(0), trid(0), pid(0),
      plane(0), module(0), 
      edep(0), time(0), xpos(0), ypos(0), 
      px(0), py(0), pz(0), 
      sizex(0), sizey(0), startx(0), starty(0)
      {}
    virtual ~GEMSimHit_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  struct HitData_t : public VDetData_t {
    UInt_t          nhits;
    std::vector<short>   *chan;           // channel (PMT)
    std::vector<unsigned int> *dataword;  // data word (encoded ADC/TDC value)
    std::vector<int>     *adc;            // ped sub ADC value
    std::vector<int>     *tdc_l;          // lead TDC (centered on t = 0)
    std::vector<int>     *tdc_t;          // trail TDC (centered on t = 0)
    bool fReadADC;                        // bool to read pedsub ADC branch
    bool fReadTDC;                        // bool to read zero centered TDC branch
    // these flags are set in constructor - just below
    HitData_t(bool readADC = true, bool readTDC = true) : nhits(0),
      chan(0), dataword(0), adc(0), tdc_l(0), tdc_t(0),
      fReadADC(readADC), fReadTDC(readTDC)
      {}
    virtual ~HitData_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  struct SampHitData_t : public HitData_t {
    //"nwords" and "adcsum" are covered with "dataword" and "adc"
    std::vector<unsigned int> *nsamps;
    std::vector< std::vector<int> > *samps_adc; // decoded ADC samples
    std::vector< std::vector<unsigned int> > *samps_datawords; // encoded ADC samples
    SampHitData_t(bool readTDC = true) : 
      HitData_t(true, readTDC),
	nsamps(0), samps_adc(0), samps_datawords(0)
      {}
    virtual ~SampHitData_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  //U for "universal"...
  struct UHitData_t : public VDetData_t {
    UInt_t          nhits;
    std::vector<short>   *chan;           // channel (PMT)
    std::vector<unsigned int> *dataword;  // data word (encoded ADC/TDC value)
    std::vector<int>     *adc;            // ped sub ADC value
    std::vector<int>     *tdc_l;          // lead TDC (centered on t = 0)
    std::vector<int>     *tdc_t;          // trail TDC (centered on t = 0)
    //variables specific to hits with samples;
    std::vector<unsigned int> *nsamps;    // actual number of samples
    std::vector< std::vector<int> > *samps_adc; // decoded ADC samples
    std::vector< std::vector<unsigned int> > *samps_datawords; // encoded ADC samples
    bool fReadADC;                        // bool to load ADC branch
    bool fReadTDC;                        // bool to load zero centered TDC branch
    bool fReadSamples;                    // bool to load samples
    // these flags are set in constructor - just below
  UHitData_t(bool readadc = true, 
	     bool readtdc = true,
	     bool readsamp = false) : nhits(0),
      chan(0), dataword(0), adc(0), tdc_l(0), tdc_t(0),
      nsamps(0), samps_adc(0), samps_datawords(0), 
      fReadADC(readadc), fReadTDC(readtdc), fReadSamples(readsamp)
      {}
    virtual ~UHitData_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
    
  const std::string kProj_str[2] = {"x", "y"};
}//end of namespace

#endif // #ifdef digsim_data_cxx
