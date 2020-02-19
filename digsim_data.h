#ifndef digsim_data_h
#define digsim_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
class TTree;


namespace SBSDigSim {

  // Purely virtual data structure for detector data in a G4SBS ROOT tree.
  struct VDetData_t {
    VDetData_t() {};
    virtual ~VDetData_t(){};
    // All sub-classes *must* implement a concrete instance of this
    virtual bool SetupBranches(TTree *t, const char* prefix) = 0;
    template<typename T>
    int SetupBranch(TTree* tree,const char* prefix,const char* varname,T &var);
  };
  
  struct PMTSimHit_t : public VDetData_t {
    UInt_t          nsimhits;
    std::vector<short>   *src;
    std::vector<short>   *trid;
    std::vector<int>     *pid;
    std::vector<short>   *chan;
    std::vector<double>  *edep;
    std::vector<int>     *npe;
    std::vector<double>  *time;
    std::vector<double>  *t_lead;
    std::vector<double>  *t_trail;
  PMTSimHit_t() : nsimhits(0),
      src(0), trid(0), pid(0), 
      chan(0), edep(0), npe(0), 
      time(0), t_lead(0), t_trail(0)
      {}
    virtual ~PMTSimHit_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  struct GEMSimHit_t : public VDetData_t {
    UInt_t          nsimhits;
    std::vector<short>   *src;
    std::vector<short>   *trid;
    std::vector<int>     *pid;
    std::vector<short>   *plane;
    std::vector<short>   *module;
    std::vector<double>  *edep;
    std::vector<double>  *time;
    std::vector<double>  *xpos;
    std::vector<double>  *ypos;
    std::vector<double>  *px;
    std::vector<double>  *py;
    std::vector<double>  *pz;
    std::vector<short>   *sizex;
    std::vector<short>   *sizey;
    std::vector<short>   *startx;
    std::vector<short>   *starty;
  GEMSimHit_t() : nsimhits(0),
      src(0), trid(0), pid(0),
      plane(0), module(0), 
      edep(0), time(0), xpos(0), ypos(0), 
      px(0), py(0), pz(0), 
      sizex(0), sizey(0), startx(0), starty(0)
      {}
    virtual ~GEMSimHit_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  struct PMTData_t : public VDetData_t {
    UInt_t          nhits;
    std::vector<short>   *chan;
    std::vector<unsigned int> *dataword;
    std::vector<int>     *adc;
    std::vector<int>     *tdc;
  PMTData_t() : nhits(0),
      chan(0), dataword(0), adc(0), tdc(0)
      {}
    virtual ~PMTData_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  struct SampPMTData_t : public PMTData_t {
    //std::vector<unsigned int> *nwords; => dataword
    //std::vector<int>     *adcsum; => adc
    std::vector< std::vector<int> > *samps_adc;
    std::vector< std::vector<unsigned int> > *samps_datawords;
  SampPMTData_t() : 
    PMTData_t(),
      samps_adc(0), samps_datawords(0)
      {}
    virtual ~SampPMTData_t(){};
    virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  
  struct GEMData_t : public VDetData_t {
    UInt_t          nhits;
    std::vector<short>   *plane;
    std::vector<short>   *module;
    std::vector<short>   *proj;
    std::vector<unsigned int> *nwords;
    std::vector< std::vector<short> > *strip;
    std::vector< std::vector<short> > *samp;
    std::vector< std::vector<int> > *samps_adc;
  GEMData_t() : nhits(0),
       plane(0), module(0), proj(0), 
       nwords(0), strip(0), samp(0), samps_adc(0) 
       {}
     virtual ~GEMData_t(){};
     virtual bool SetupBranches(TTree *t, const char* prefix);
  };
  
  
  
}

#endif // #ifdef digsim_data_cxx
