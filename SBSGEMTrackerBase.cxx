#include "SBSGEMTrackerBase.h"
#include "SBSGEMModule.h"
#include "TRotation.h"

SBSGEMTrackerBase::SBSGEMTrackerBase(){ //Set default values of important parameters: 
  Clear();

  fIsMC = false;
  fNmodules = 0;
  fNlayers = 0;
  fTrackingAlgorithmFlag = 2;

  fMinHitsOnTrack = 3;

  fOnlinePedestalSubtraction = true;
  fZeroSuppress = true;
  fZeroSuppressRMS = 10.0; //10 ADC channels. We are free to define how this is actually used;

  fGridBinWidthX = 0.01; //1 cm = 10 mm;
  fGridBinWidthY = 0.01; //1 cm = 10 mm;

  fTrackChi2Cut = 100.0; //Max. chi2/ndf for a combination of hits to form a track

  // set defaults for constraint points and constraint widths:
  fConstraintPoint_Front.SetXYZ(0,0,0);
  fConstraintPoint_Back.SetXYZ(0,0,10.0);

  //wide-open constraints for now:
  fConstraintWidth_Front.Set( 1.5, 0.5 );
  fConstraintWidth_Back.Set( 1.5, 0.5 ); 
}

void SBSGEMTrackerBase::~SBSGEMTrackerBase(){
  //for now, do nothing; let the derived classes handle the clearing out of the modules
  
}

void SBSGEMTrackerBase::Clear(){ //Clear out any event-specific stuff
  fNtracks_found = 0;
  fNhitsOnTrack.clear();
  fModListTrack.clear();
  fHitListTrack.clear();
  fresidu_hits.clear();
  fresidv_hits.clear();
  feresidu_hits.clear();
  feresidv_hits.clear();

  fXtrack.clear();
  fYtrack.clear();
  fXptrack.clear();
  fYptrack.clear();
  fChi2Track.clear();

  fclustering_done = false;
  ftracking_done = false;
  
}

void SBSGEMTrackerBase::InitLayerCombos() { //It is assumed that this will be called by the ReadDatabase/Init method of the derived classes after loading the necessary information from the database:

  //Just in case:
  fLayerCombinations.clear();
  
  //Initialize the list of layer combinations by total number of layers on the combo:
  for( int icombo=0; icombo<pow(2,fNlayers), icombo++ ){ //loop over all possible combinations of layers:
    
    vector<int> layercombo; //temporary array to hold list of layers fired in combo
    int nlayersoncombo=0; //count number of fired layers in combo
    for( int ilayer=0; ilayer<fNlayers; ilayer++ ){ //loop over all layers
      int testbit = pow(2,ilayer); 
      if( testbit & icombo != 0 ){ //bitwise AND of icombo and 2^ilayer nonzero:
	nlayersoncombo++; //this layer is on the combo
	layercombo.push_back( ilayer ); //add it to the list on this combo
      }
    }

    if( nlayersoncombo >= fMinHitsOnTrack ){ //Only consider combinations of layers with at least the minimum number required
      fLayerCombinations[nlayersoncombo].push_back( layercombo ); //Add this combo to the list of combos, mapped by the number of layers in the combination:
    }
    
  }
}

void SBSGEMTrackerBase::InitGridBins() {
  //we assume that the database has already been read and the module geometry is already specified here:
  //Loop over layers and modules within each layer, and set the size of the active area by layer:

  //clear out any existing data, just in case:
  fXmin_layer.clear();
  fXmax_layer.clear();
  fYmin_layer.clear();
  fYmax_layer.clear();

  fGridNbinsX_layer.clear();
  fGridNbinsY_layer.clear();
  //
  fGridXmin_layer.clear();
  fGridXmax_layer.clear();
  fGridYmin_layer.clear();
  fGridYmax_layer.clear();
  
  for( int layer = 0; layer<fNlayers; layer++ ){
    std::set<int> modlist_layer = fModuleListByLayer[layer];

    //Initialize grid active area for this layer:
    double xgmin_all = 1.e12, ygmin_all = 1.e12, xgmax_all = -1.e12, ygmax_all = -1.e12; 

    double zsum = 0.0;
    
    for( auto imod = modlist_layer.begin(); imod != modlist_layer.end(); ++imod ){
      
      int module = *imod; //if the database is sensibly constructed, this should refer to the index in the array of modules:

      SBSGEMModule *mtemp = fModules[module];

      //Get origin coordinates:
      TVector3 modpos = mtemp->GetOrigin();

      //computation of average Z coordinate of modules in this layer:
      zsum += modpos.Z();
      
      //Get half-width of module along X and Y:
      double Lx_mod = mtemp->GetXSize()/2.0;
      double Ly_mod = mtemp->GetYSize()/2.0;

      //get positions of the four corners of the active area (which is assumed rectangular for SBS GEMs):
      TVector3 Corner1 = modpos - Lx_mod * mtemp->GetXax() - Ly_mod * mtemp->GetYax();
      TVector3 Corner2 = modpos + Lx_mod * mtemp->GetXax() - Ly_mod * mtemp->GetYax();
      TVector3 Corner3 = modpos - Lx_mod * mtemp->GetXax() + Ly_mod * mtemp->GetYax();
      TVector3 Corner4 = modpos + Lx_mod * mtemp->GetXax() + Ly_mod * mtemp->GetYax();

      //Check all four corners even though ONLY corners 1 and 3 are likely to define the minimum X
      xgmin_all = Corner1.X() < xgmin_all ? Corner1.X() : xgmin_all;
      xgmin_all = Corner2.X() < xgmin_all ? Corner2.X() : xgmin_all;
      xgmin_all = Corner3.X() < xgmin_all ? Corner3.X() : xgmin_all;
      xgmin_all = Corner4.X() < xgmin_all ? Corner4.X() : xgmin_all;

      //Check all four corners even though ONLY corners 2 and 4 are likely to define the maximum X
      xgmax_all = Corner1.X() > xgmax_all ? Corner1.X() : xgmax_all;
      xgmax_all = Corner2.X() > xgmax_all ? Corner2.X() : xgmax_all;
      xgmax_all = Corner3.X() > xgmax_all ? Corner3.X() : xgmax_all;
      xgmax_all = Corner4.X() > xgmax_all ? Corner4.X() : xgmax_all;

      //Check all four corners even though ONLY corners 1 and 2 are likely to define the minimum Y
      ygmin_all = Corner1.Y() < ygmin_all ? Corner1.Y() : ygmin_all;
      ygmin_all = Corner2.Y() < ygmin_all ? Corner2.Y() : ygmin_all;
      ygmin_all = Corner3.Y() < ygmin_all ? Corner3.Y() : ygmin_all;
      ygmin_all = Corner4.Y() < ygmin_all ? Corner4.Y() : ygmin_all;

      //Check all four corners even though ONLY corners 3 and 4 are likely to define the maximum Y
      ygmax_all = Corner1.Y() > ygmax_all ? Corner1.Y() : ygmax_all;
      ygmax_all = Corner2.Y() > ygmax_all ? Corner2.Y() : ygmax_all;
      ygmax_all = Corner3.Y() > ygmax_all ? Corner3.Y() : ygmax_all;
      ygmax_all = Corner4.Y() > ygmax_all ? Corner4.Y() : ygmax_all;

      fXmin_layer[layer] = xgmin_all;
      fXmax_layer[layer] = xgmax_all;
      fYmin_layer[layer] = ygmin_all;
      fYmax_layer[layer] = ygmax_all;

    } //end loop over list of modules in this layer

    fZavgLayer[layer] = zsum/double( modlist_layer.size() );

    fGridXmin_layer[layer] = fXmin_layer[layer] - 0.5*fGridBinWidthX;
    int nbinsx = 0;
    while( fGridXmin_layer[layer] + nbinsx * fGridBinWidthX < fXmax_layer[layer] + 0.5*fGridBinWidthX ){
      nbinsx++;
    }

    fGridXmax_layer[layer] = fGridXmin_layer[layer] + nbinsx * fGridBinWidthX;
    fGridNbinsX_layer[layer] = nbinsx;

    fGridYmin_layer[layer] = fYmin_layer[layer] - 0.5*fGridBinWidthY;
    int nbinsy = 0;
    while( fGridYmin_layer[layer] + nbinsy * fGridBinWidthY < fYmax_layer[layer] + 0.5*fGridBinWidthY ){
      nbinsy++;
    }

    fGridYmax_layer[layer] = fGridYmin_layer[layer] + nbinsy * fGridBinWidthY;
    fGridNbinsY_layer[layer] = nbinsy;
    
  } //end loop over layers
  
}

void SBSGEMTrackerBase::hit_reconstruction(){

  fclustering_done = true;
  
  //Loop over all the GEM modules and invoke their cluster-finding methods with search region constraints:
  for( int imodule=0; imodule<fNmodules; imodule++ ){
    SBSGEMModule *mod = fModules[imodule];

    if( !fUseConstraint ){ //call find_2D hits for the module without any constraint:
      mod->find_2Dhits();
    } else { //Determine constraint center and constraint width at this module:
      TVector3 modpos = mod->GetOrigin();
      TVector3 modzaxis = mod->GetZax();

      TVector3 constraint_direction = (fConstraintPoint_Back - fConstraintPoint_Front).Unit();

      //Calculate intersection point with plane of module:
      // recall modzaxis dot ( r - modpos ) = 0
      // r = r0 + s * nhat
      // modzaxis dot ( r0 + s* nhat - modpos ) = 0
      // s (modzaxis dot nhat) = modzaxis dot (modpos - r0) --> s = modzaxis dot (modpos - r0)/modzaxis dot nhat:

      //This is the distance along the line defined by the front and rear constraint points to its intersection point with the plane of the module:
      double sintersect = modzaxis.Dot( modpos - fConstraintPoint_Front )/( modzaxis.Dot( constraint_direction ) );

      //The following is the intersection point of the line defined by the front and rear constraint points with the plane of the module:
      TVector3 constraint_intersect = fConstraintPoint_Front + sintersect * constraint_direction;

      //Define the constraint width at each module via linear interpolation between the front and rear constraint widths:
      //Note: It is assumed that the front and rear constraint points and widths will be defined to be just outside the entire physical z extent
      //      of all the layers,
      //      such that every module lies between the front and rear constraint points (along z) by definition, and therefore that
      //      "interp_frac" below lies between zero and one in virtually all cases
      //      
      
      double interp_frac = sintersect / (fConstraintPoint_Back-fConstraintPoint_Front).Mag();
      
      TVector2 constraint_width_module( fConstraintWidth_Front.X() * (1.-interp_frac) + fConstraintWidth_Back.X() * interp_frac,
					fConstraintWidth_Front.Y() * (1.-interp_frac) + fConstraintWidth_Back.Y() * interp_frac );

      //compute constraint in "local" module coordinates:
      TVector3 constraint_intersect_local = TrackToDetCoord( constraint_intersect );
      
      TVector2 constraint_center_module( constraint_intersect_local.X(), constraint_intersect_local.Y() );

      //Do 2D hit reconstruction in the search region defined by the constraints at this module.
      //First check if any part of the search region overlaps the module active area:
      
      mod->find_2Dhits( constraint_center_module, constraint_width_module );
    }
      
  }
}

// Standard, "default" fast track-finding algorithm (based on SBSGEM_standalone code by Andrew Puckett):
void SBSGEMTrackerBase::find_tracks(){ 

  fNtracks_found = 0;
  
  if( !fclustering_done ){ //This shouldn't be called before hit reconstruction, but if it is, then this routine can call the hit reconstruction:
    hit_reconstruction();
  }
  
  ftracking_done = true;
  //It is assumed that when we reach this stage, the hit reconstruction will have already been called. 
  
  //Define and initialize the "hit lists" and other local arrays that we will need to do the tracking:
  //These "local hit list" arrays will NOT be modified throughout the track-finding procedure:
  std::set<int> layers_with_2Dhits; //list of tracking layers with at least one 2D hit reconstructed (within the track search region)
  std::map<int,int> N2Dhits_layer; // key = layer, mapped value = total number of reconstructed 2D hits
  std::map<int,std::vector<int> > modindexhit2D; //key = layer, mapped value = module index of hits in that layer:
  std::map<int,std::vector<int> > clustindexhit2D; //key = layer, mapped value = index of hits within 2D cluster array of module in question
  std::map<int,std::vector<bool> > hitused2D; //flag to tell each track-finding iteration whether hit was already used in a previous track

  for( int imodule=0; imodule<fModules.size(); imodule++ ){ //loop over all the 2D hits in all modules (track search region was already enforced in hit_reconstruction)
    int layer = fModules[imodule]->GetLayer();

    int n2Dhits_mod = fModules[imodule]->fHits.size(); 
    
    for( int ihit=0; ihit<n2Dhits; ihit++ ){
      sbsgemhit_t hittemp = fModules[imodule]->fHits[ihit];

      if( hittemp.keep ){
	layers_with_2Dhits.insert( layer );
	modindexhit2D[layer].push_back( imodule ); //module index
	clustindexhit2D[layer].push_back( ihit ); //index in this module's cluster array
	hitused2D[layer].push_back( false );

	N2Dhits_layer[layer] = modindexhit2D[layer].size();
      }
    }
  }

  //At this stage the local "hit lists" that we need for the tracking are initialized. Let's get started:

  if( layers_with_2Dhits.size() >= fMinHitsOnTrack ){ //Then we have enough layers to do tracking:
    bool foundtrack = true;
    while( foundtrack ){ //As long as we find a track on each iteration, keep looking for more tracks:
      foundtrack = false;

      int nhitsrequired = layers_with_2Dhits.size(); //initially we favor tracks with the largest possible number of hits; if we fail to find a track at this hit requirement, we decrement the number of required hits as long as it exceeds the minimum

      while( nhitsrequired >= fMinHinHitsOnTrack ){ //as long as the current hit requirement exceeds the minimum, 
	foundtrack = false;

	
	
	if( !foundtrack ){
	  nhitsrequired--;
	}
      }
      
    } //end while( foundtrack )


  } //end check of sufficient layers with hits to do tracking

  
  

}
