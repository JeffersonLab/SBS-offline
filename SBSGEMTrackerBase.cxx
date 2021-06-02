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

  fOnlineZeroSuppression = false;
  fZeroSuppress = true;
  fZeroSuppressRMS = 5.0; //5 sigma

  fGridBinWidthX = 0.01; //1 cm = 10 mm;
  fGridBinWidthY = 0.01; //1 cm = 10 mm;

  fGridEdgeToleranceX = 0.003; //3 mm default
  fGridEdgeToleranceY = 0.003; //3 mm default
  
  fTrackChi2Cut = 100.0; //Max. chi2/ndf for a combination of hits to form a track

  fSigma_hitpos = 0.0001; //100 um
  fSigma_hitshape = 0.0004; //0.4 mm
  
  // set defaults for constraint points and constraint widths:
  fConstraintPoint_Front.SetXYZ(0,0,-10.0);
  fConstraintPoint_Back.SetXYZ(0,0,10.0);

  //wide-open constraints for now (the units are meters here (mm or cm would be more natural but whatever))
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

//This is called by the Init() methods of derived classes:
void SBSGEMTrackerBase::CompleteInitialization(){
  fLayers.clear();
  fLayerByIndex.clear();
  fNumModulesByLayer.clear();
  fModuleListByLayer.clear();
  //loop on the array of (initialized) modules and fill out missing info:
  
  for( int imod=0; imod<fModules.size(); imod++ ){
    int layer = fModules[imod]->fLayer;

    auto newlayer = fLayers.insert( layer );

    fModuleListByLayer[layer].insert( imod );
  }

  fNmodules = fModules.size();
  fNlayers = fLayers.size();

  for( auto ilay = fLayers.begin(); ilay != fLayers.end(); ++ilay ){
    int layer = *ilay;
    fLayerByIndex.push_back( layer ); //this is a vector version of the layer list, unless the user has defined something weird, layer[layerindex] = layerindex for layerindex = 0, ..., Nlayers-1
    fNumModulesByLayer[layer] = fModuleListByLayer[layer].size();
  }

  //make sure the user has defined something sensible:
  fMinHitsOnTrack = std::max(3,std::min(fNlayers,fMinHitsOnTrack) );
  
  InitLayerCombos();
  InitGridBins();
}

//This only needs to be done ONCE (after loading the geometry from the database!)
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
  
  for( int ilayer = 0; ilayer<fNlayers; ilayer++ ){
    int layer = fLayerByIndex[ilayer];
    
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

//Initialize the "hit list" arrays that are used by the track-finding algorithm: these arrays are UNCHANGING throughout the iterations of track-finding:
void SBSGEMTrackerBase::InitHitList(){

  //clear out any old information:
  layers_with_2Dhits.clear();
  N2Dhits_layer.clear();
  modindexhit2D.clear();
  clustindexhit2D.clear();
  hitused2D.clear();
  
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
}

void SBSGEMTrackerBase::InitFreeHitList(){
  //We should clear these things out at the beginning of each iteration just in case:
  layerswithfreehits.clear();
  freehitlist_layer.clear();
  freehitcounter.clear();
  Nfreehits_layer.clear();
  
  Nfreehits_binxy_layer.clear();
  freehitlist_binxy_layer.clear();
  
  //This is called at the beginning of each track-finding iteration:
  for( auto ilay = layers_with_2Dhits.begin(); ilay != layers_with_2Dhits.end(); ++ilay ){
    int layer = *ilay;
    Nfreehits_layer[layer] = 0;
    
    int nbins_gridxy = fGridNbinsX_layer[layer]*fGridNbinsY_layer[layer];

    //resize the free hit counter and free hit list for the 2D grid bins to the number of grid bins:
    Nfreehits_binxy_layer[layer].resize( nbins_gridxy );
    freehitlist_binxy_layer[layer].resize( nbins_gridxy );

    //loop over the 2D grid bins and initialize free hit counter and freehit list:
    for( int bin=0; bin<nbins_gridxy; bin++ ){
      Nfreehits_binxy_layer[layer][bin] = 0;
      freehitlist_binxy_layer[layer][bin].clear();
    }

    //Now loop over all the hits and fill up the free hit list arrays:
    for( int ihit=0; ihit<N2Dhits_layer[layer]; ihit++ ){
      //Make sure this hit is not already used in a track:
      if( !hitused2D[layer][ihit] ){ //
	Nfreehits_layer[layer]++;
	freehitlist_layer[layer].push_back( ihit ); //recall that "ihit" is the index in the unchanging "hit list" array

	int module = modindexhit2D[layer][ihit]; 
	int clustidx = clustindexhit2D[layer][ihit]; //index in module hit array
	
	//Next, let's grab the global X and Y coordinates of the hit and assign it to a grid bin:
	double xgtemp = fModules[module]->fHits[clustidx].xghit;
	double ygtemp = fModules[module]->fHits[clustidx].yghit;

	//Now figure out which X and Y bin of the 2D grid this hit belongs in:

	int binxtemp = int( (xgtemp - fGridXmin_layer[layer])/fGridBinWidthX );
	int binytemp = int( (ygtemp - fGridYmin_layer[layer])/fGridBinWidthY );

	//check that hit coordinates are within the limits of the grid (in principle this check is unnecessary since the grid is defined
	// to contain the entire active area of the layer, but just to avoid any potential seg. faults we check anyway):
	if( binxtemp >= 0 && binxtemp < fGridNbinsX_layer[layer] &&
	    binytemp >= 0 && binytemp < fGridNbinsY_layer[layer] ){
	  int binxytemp = binxtemp + fGridNbinsX_layer[layer]*binytemp;

	  Nfreehits_binxy_layer[layer][binxytemp]++;
	  freehitlist_binxy_layer[layer][binxytemp].push_back( ihit ); //Here again, ihit locates this hit within the unchanging "hit list" array
	  
	}
	
      } //check hit not already used in track
    } //loop over all hits in layer

    if( Nfreehits_layer[layer] > 0 ){
      layerswithfreehits.insert(layer);
      freehitcounter[layer] = 0;
    }
    
  } //loop over all layers
}

void SBSGEMTrackerBase::hit_reconstruction(){

  fclustering_done = true;
  
  //Loop over all the GEM modules and invoke their cluster-finding methods with search region constraints:
  for( int imodule=0; imodule<fNmodules; imodule++ ){
    SBSGEMModule *mod = fModules[imodule];

    if( !fUseConstraint ){ //call find_2D hits for the module without any constraint:
      mod->find_2Dhits();
    } else {
      
      TVector3 constraint_direction = (fConstraintPoint_Back - fConstraintPoint_Front).Unit();

      //Calculate intersection point with plane of module:
      // recall modzaxis dot ( r - modpos ) = 0
      // r = r0 + s * nhat
      // modzaxis dot ( r0 + s* nhat - modpos ) = 0
      // s (modzaxis dot nhat) = modzaxis dot (modpos - r0) --> s = modzaxis dot (modpos - r0)/modzaxis dot nhat

      //The following is the intersection point of the line defined by the front and rear constraint points with the plane of the module:
      double sintersect; //dummy variable to hold distance along track to intersection from front constraint point:
      TVector3 constraint_intersect = TrackIntersect( imodule, fConstraintPoint_Front, constraint_direction, sintersect );
      
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
      TVector3 constraint_intersect_local = mod->TrackToDetCoord( constraint_intersect );
      
      TVector2 constraint_center_module( constraint_intersect_local.X(), constraint_intersect_local.Y() );

      //Do 2D hit reconstruction in the search region defined by the constraints at this module.
      //First check if any part of the search region overlaps the module active area:
      
      mod->find_2Dhits( constraint_center_module, constraint_width_module );
    }
      
  }
}

// Standard "fast" track-finding algorithm (based on SBSGEM_standalone code by Andrew Puckett):
void SBSGEMTrackerBase::find_tracks(){ 

  //should this method invoke clear()? Yes: Clear() just clears out all the track arrays. It is assumed that this method will only be called once per event.
  //Although that is probably not correct; it might be called as many as two times. Anyway, for now, let's use it, might need to revisit later:
  //Clear();
  
  fNtracks_found = 0;
  
  if( !fclustering_done ){ //This shouldn't be called before hit reconstruction, but if it is, then this routine will call the hit reconstruction:
    hit_reconstruction();
  }
  
  ftracking_done = true;
  //It is assumed that when we reach this stage, the hit reconstruction will have already been called. 

  //Initialize the (unchanging) hit list that will be used by the rest of the tracking procedure:
  InitHitList();

  //At this stage the static "hit lists" that we need for the tracking are initialized. Let's get started:

  if( layers_with_2Dhits.size() >= fMinHitsOnTrack ){ //Then we have enough layers to do tracking:
    //bool foundtrack = true; rendered unnecessary by the removal of the outermost, redundant while loop:
      
    int nhitsrequired = layers_with_2Dhits.size(); //initially we favor tracks with the largest possible number of hits; if we fail to find a track at this hit requirement, we decrement the number of required hits as long as it exceeds the minimum

    while( nhitsrequired >= fMinHinHitsOnTrack ){ //as long as the current minimum hit requirement exceeds the minimum hits to define a track, we look for more tracks with
      // nhitsrequired hits:
      bool foundtrack = false;

      //The goal will be to loop over all possible combinations of layers at the current hit requirement, and find the combination of one hit per layer with minimum chi2:
      //At the beginning of each iteration of this loop, we need to populate several more arrays, the "free hit list" mapped by layer, and
      //also the "free hit list" mapped by layer and 2D grid bin:

      // We moved the initialization of the "free hit" list for each track-finding iteration to its own function,
      // and we moved the relevant arrays to data members of the class

      //This happens once per track-finding iteration: if any tracks were found on the previous iteration, then their hits (and any others sharing the same 1D U or V clusters)
      //will have been marked as used, reducing the number of "available" hits for finding additional tracks:
      InitFreeHitList(); 

      if( layerswithfreehits.size() >= nhitsrequired ){ //check that the number of layers with free hits is at least equal to the current minimum hit requirement:
	//The basic algorithm should do the following:

	// 1. Loop over all possible layer combinations for which the number of layers is equal to the current minimum hit requirement:
	// 2. Within each layer combination at the current minimum hit requirement, check that every layer has at least one free hit. If so, proceed:
	// 3. For all possible combinations of one hit from each of the two outermost layers, calculate the straight line between the two points, and project the
	//    straight line to each layer
	// 4. At each intermediate layer, populate the list of hits in the grid bins pointed to by the candidate track. If the projected track is close to the edge of
	//    the bin in either X or Y, also consider hits in the adjacent X and/or Y bins
	// 5. IFF every intermediate layer contains at least one hit in the grid bin(s) sufficiently close to the projected track, proceed to loop over the
	//    hits in the intermediate layers, using either the "fast" method (find the hit in each layer closest to the projected track) or the "brute force" method
	//    (loop on all possible hit combinations in the intermediate layers using "odometer" algorithm), which is probably more accurate, but slower
	// 6. For each candidate hit combination, calculate the chi2 of a straight line fit.

	// The actual loop over hit combinations starts here, define local variables needed to store best hit combination, best track and residuals, and minimum chi2:
	bool firstgoodcombo = true;
	map<int,int> besthitcombo;
	double minchi2 = 1.e20; //arbitrary large number initially

	vector<double> besttrack(4); //x, y, x', y'

	//Let's carry around the track fitting residuals so that we don't have to repeat the calculation when adding the fitted track with best chi2 to the track arrays:
	vector<double> uresidbest, vresidbest;
	  
	for( int icombo=0; icombo<fLayerCombinations[nhitsrequired].size(); icombo++ ){
	  int minlayer = fNlayers + 1;
	  int maxlayer = -1;

	  //list of layers to test on this hit combination (all layers have to fire in order to proceed):
	  set<int> layerstotest;

	  //Also record outermost layers for fast track-finding using "grid search":
	    
	  for( int ihit=0; ihit<nhitsrequired; ihit++ ){
	    int layeri = fLayerCombinations[nhitsrequired][icombo][ihit];
	    int layer = fLayerByIndex[layeri];
	    if( layerswithfreehits.find( layer ) != layerswithfreehits.end() ){ //check that this layer has unused hits:
	      layerstotest.insert( layer );

	      minlayer = (layer < minlayer) ? layer : minlayer;
	      maxlayer = (layer > maxlayer) ? layer : maxlayer;
	    }
	  }

	  //Calculate the total number of possible combinations of one hit from each of the two outermost layers:
	  long ncombos_minmax = Nfreehits_layer[minlayer]*Nfreehits_layer[maxlayer];
	    
	  if( layerstotest.size() < nhitsrequired ){
	    //skip this layer combination if any layers lack free hits at the current minimum hit requirement:
	    continue;
	  }

	  // If the number of hit combinations in the two outermost layers exceeds the maximum, try to find an alternate pair of layers
	  // with hit combinations below the maximum:
	  if( ncombos_minmax > fMaxHitCombinations ){
	    // if the two outermost layers in this combo have too many hit combinations, find the combination of layers
	    // with the largest lever arm in z such that the number of combinations is less than the
	    // maximum:
	    int maxdiff = 0;
	      
	    for( auto ilay = layerstotest.begin(); ilay != layerstotest.end(); ++ilay ){
	      int layeri = *ilay;
	      for( auto jlay = ilay; jlay != layerstotest.end(); ++jlay ){
		int layerj = *jlay;
		if( layerj > layeri ){
		  long ncombostest = Nfreehits_layer[layeri]*Nfreehits_layer[layerj];

		  if( ncombostest >= 1 && ncombostest <= fMaxHitCombinations &&
		      layerj - layeri > maxdiff ){
		    minlayer = layeri;
		    maxlayer = layerj;
		    maxdiff = layerj - layeri;
		    ncombos_minmax = ncombostest;
		  }	  
		}
	      }	
	    }
	  } //end check of ncombos_minmax < maxhitcombinations

	    //If the number of hit combinations STILL exceeds the maximum, skip this layer combination:
	  if( ncombos_minmax > fMaxHitCombinations ){
	    continue;
	  }
	    
	  //Next: loop over all hits in the two outermost layers and form track from each combination:
	  for( int ihit = 0; ihit<Nfreehits_layer[minlayer]; ihit++ ){
	    for( int jhit = 0; jhit<Nfreehits_layer[maxlayer]; jhit++ ){
	      //The track search region constraint should have already been enforced at the 2D hit reconstruction stage, so additional checks here are probably unnecessary.
		
	      int hitmin = freehitlist_layer[minlayer][ihit];
	      int hitmax = freehitlist_layer[maxlayer][jhit];

	      int modmin = modindexhit2D[minlayer][hitmin];
	      int modmax = modindexhit2D[maxlayer][hitmax];

	      int clustmin = clustindexhit2D[minlayer][hitmin];
	      int clustmax = clustindexhit2D[maxlayer][hitmax];

	      //Get 3D global coordinates of the two hits:
	      TVector3 hitpos_min = GetHitPosGlobal( modmin, clustmin );
	      TVector3 hitpos_max = GetHitPosGlobal( modmax, clustmax );

	      // populate the list of layers other than minlayer and maxlayer to build the track:
	      std::set<int> otherlayers;

	      for( auto ilay = layerstotest.begin(); ilay != layerstotest.end(); ++ilay ){
		int thislayer = *ilay;
		if( thislayer != minlayer && thislayer != maxlayer ){
		  otherlayers.insert( *ilay );
		}
	      }

	      //This array will hold the list of free hits in layers other than minlayer and maxlayer falling in 2D grid bins
	      //close to the track projection:
	      std::map<int,std::vector<int> > freehitlist_otherlayers_goodxy; 
		
	      //The next step is to calculate the straight line passing through the two points from minlayer and maxlayer:
	      // double xptrtemp = (hitpos_max.X() - hitpos_min.X())/(hitpos_max.Z()-hitpos_min.Z());
	      // double yptrtemp = (hitpos_max.Y() - hitpos_min.Y())/(hitpos_max.Z()-hitpos_min.Z());
		
		
	      //Project track to z = 0 plane:
	      double xptrtemp = (hitpos_max.X() - hitpos_min.X())/(hitpos_max.Z() - hitpos_min.Z() );
	      double yptrtemp = (hitpos_max.Y() - hitpos_min.Y())/(hitpos_max.Z() - hitpos_min.Z() );

	      //Track coordinates at Z = 0:
	      double xtrtemp = 0.5*( hitpos_max.X() - xptrtemp * hitpos_max.Z() + hitpos_min.X() - xptrtemp * hitpos_min.Z() );
	      double ytrtemp = 0.5*( hitpos_max.Y() - yptrtemp * hitpos_max.Z() + hitpos_min.Y() - yptrtemp * hitpos_min.Z() );

	      //Next we will project the track to the average Z coordinate of each layer in "otherlayers" and check for hits in nearby grid bins:

	      bool nextcomboexists = true;

	      //clear out the "free hit" counter for looping over combinations:
	      freehitcounter.clear();
		
	      for( auto ilay = otherlayers.begin(); ilay != otherlayers.end(); ++ilay ){
		int layer = *ilay;

		double xproj = xtrtemp + xptrtemp * fZavgLayer[layer];
		double yproj = ytrtemp + yptrtemp * fZavgLayer[layer];

		int binxtemp = int( (xproj - fGridXmin_layer[layer])/fGridBinWidthX );
		int binytemp = int( (yproj - fGridYmin_layer[layer])/fGridBinWidthY );

		//check if this bin is in-range:

		if( binxtemp >= 0 && binxtemp < fGridNbinsX_layer[layer] &&
		    binytemp >= 0 && binytemp < fGridNbinsY_layer[layer] ){

		  int binxhi = binxtemp;
		  int binxlo = binxtemp;
		  int binyhi = binytemp;
		  int binylo = binytemp;

		  //the following two variables are the differences between the track projection to this layer and the low edge of the bin in both X and Y:
		  double binxdiff = xproj - (fGridXmin_layer[layer]+binxtemp*fGridBinWidthX);
		  double binydiff = yproj - (fGridYmin_layer[layer]+binytemp*fGridBinWidthY);
		    
		  //If x or y projection is close to the low edge of bin, include the neighboring bin on the low side in the analysis, assuming it exists:
		  if( binxdiff < fGridEdgeToleranceX && binxtemp > 0 ) binxlo = binxtemp-1;
		  if( binydiff < fGridEdgeToleranceY && binytemp > 0 ) binylo = binytemp-1;
		  //I x or y projection is close to the high edge of the bin, include the neighboring bin on high side in the analysis, assuming it exists:
		  if( fGridBinWidthX - binxdiff < fGridEdgeToleranceX && binxtemp + 1 < fGridNbinsX_layer[layer] ) binxhi = binxtemp+1;
		  if( fGridBinWidthY - binydiff < fGridEdgeToleranceY && binytemp + 1 < fGridNbinsY_layer[layer] ) binyhi = binytemp+1;

		  //now loop over the relevant grid bins (up to 2 in X and Y) in this layer and fill the "reduced" free hit list:
		  for( int binx = binxlo; binx <= binxhi; binx++ ){
		    for( int biny = binylo; biny <= binyhi; biny++ ){
		      int binxy = binx + fGridNbinsX_layer[layer]*biny;

		      for( int khit=0; khit<freehitlist_binxy_layer[layer][binxy].size(); khit++ ){
			freehitlist_otherlayers_goodxy[layer].push_back( freehitlist_binxy_layer[layer][binxy][khit] );
		      }
			
		    }
		  }
		} //end check on grid bin of projected track in range

		  //This check enforces that all layers other than minlayer and maxlayer have at least one hit in the relevant 2D grid bins:
		if( freehitlist_otherlayers_goodxy.find(layer) == freehitlist_otherlayers_goodxy.end() ) nextcomboexists = false;

		freehitcounter[layer] = 0;
		  
	      } //end loop on layers other than minlayer and maxlayer

		//Next, we will loop on all possible combinations of one hit from each of the layers other than minlayer and maxlayer:
	      if( nextcomboexists ){
		bool firstcombo = true;

		std::map<int,int> hitcombo;
		  
		while( nextcomboexists = GetNextCombo( otherlayers, freehitlist_otherlayers_goodxy, freehitcounter, hitcombo, firstcombo ) ){
		  // I think that the assignment of the result of GetNextCombo() to nextcomboexists in the while loop condition renders an extra check of the value of
		  // nextcomboexists unnecessary
		  //Then we form the track from minhit, maxhit, and hitcombo, and check if this hit combination has better chi2 than any previous one:

		  //First, add the hits from minlayer and maxlayer to the combo:
		  hitcombo[minlayer] = hitmin;
		  hitcombo[maxlayer] = hitmax;

		  double xptrtemp, yptrtemp, xtrtemp, ytrtemp, chi2ndftemp;

		  //This declaration might shadow another one up above
		  //(actually it DOESN'T: the ones above are for the "best" hit combo, these are temporary dummy variables. Proceed)
		  vector<double> uresidtemp, vresidtemp;
		    
		  //Fit a track to the current hit combination:
		  //NOTE: the FitTrack method computes the line of best fit and chi2 and gives us the hit residuals:
		  FitTrack( hitcombo, xtrtemp, ytrtemp, xptrtemp, yptrtemp, chi2ndftemp, uresidtemp, vresidtemp );

		  if( firstgoodcombo || chi2ndftemp < minchi2 ){
		    firstgoodcombo = false;
		    minchi2 = chi2ndftemp;

		    besthitcombo = hitcombo;

		    //record track properties so we don't need to re-fit later
		    besttrack[0] = xtrtemp;
		    besttrack[1] = ytrtemp;
		    besttrack[2] = xptrtemp;
		    besttrack[3] = yptrtemp;

		    //Perhaps this is an inefficent copy of vector<double>, but probably fine compared to repeating the calculation of residuals later on:
		    uresidbest = uresidtemp;
		    vresidbest = vresidtemp; 
		  }

		  //clear hitcombo just so we start fresh each iteration of the loop: this is PROBABLY unnecessary, but safer than not doing so:
		  hitcombo.clear();
		    
		} //end while( nextcomboexists )
	      } //end if( nextcomboexists )		       
		
	    } //end loop over hits in maxlayer
	  } //end loop over hits in minlayer
	} //end loop over layer combinations at current minimum hit requirement

	  //We treat all layer combinations at the same minimum hit requirement on an equal footing as far as track-finding is concerned:
	if( !firstgoodcombo && minchi2 < fTrackChi2Cut ){ //then we found at least one candidate track:
	  foundtrack = true;

	  // "AddTrack" takes care of incrementing fNtracks_found
	    
	  AddTrack( besthitcombo, besttrack, minchi2, uresidbest, vresidbest );
	    
	}
	  
      } //end check on layers with free hits >= nhitsrequired
		
      if( !foundtrack ){ //If we didn't find any tracks at the current minimum hit requirement, then reduce the minimum, and see if we can find a good track
	// (or tracks) with one fewer hit. Otherwise, we search again at the current minimum hit requirement:
	nhitsrequired--;
      }
    } //end while(nhitsrequired >= minhits ) 
  } //end check of sufficient layers with hits to do tracking
}

//The next function determines the line of best fit through a combination of hits, without calculating residuals or chi2.
// Note that these equations assume all hits are to be given equal weights. You will need a different function if you want to use different weights for different hits:
void SBSGEMTrackerBase::CalcLineOfBestFit( const std::map<int,int> &hitcombo, double &xtrack, double &ytrack, double &xptrack, double &yptrack ){
  double sumx = 0.0, sumy = 0.0, sumz = 0.0, sumxz = 0.0, sumyz = 0.0, sumz2 = 0.0;
  
  int nhits = 0;
  
  for( auto ilayer=hitcombo.begin(); ilayer != hitcombo.end(); ++ilayer ){
    int layer = ilayer->first;   //layer 
    int hitidx = ilayer->second; //index in the "hit list" array
    
    //grab hit coordinates:
    int module = modindexhit2D[layer][hitidx];
    int clustidx = clustindexhit2D[layer][hitidx];
    
    TVector3 hitpos_global = GetHitPosGlobal( module, clustidx );
    
    //we don't use the u and v coordinates until the chi2 calculation, which comes later:
    //double uhit = fModules[module]->fHits[clustidx].uhit;
    //double vhit = fModules[module]->fHits[clustidx].vhit; 
    
    sumx += hitpos_global.X();
    sumy += hitpos_global.Y();
    sumz += hitpos_global.Z();
    sumxz += hitpos_global.X()*hitpos_global.Z();
    sumyz += hitpos_global.Y()*hitpos_global.Z();
    sumz2 += pow(hitpos_global.Z(),2);
    
    nhits++;
  }
  
  //now compute line of best fit:
  double denom = (sumz2 * nhits - pow(sumz,2) );
  
  xptrack = (nhits*sumxz - sumx*sumz)/denom;
  yptrack = (nhits*sumyz - sumy*sumz)/denom;
  xtrack = (sumx * sumz2 - sumxz * sumz)/denom;
  ytrack = (sumy * sumz2 - sumyz * sumz)/denom;
}

void SBSGEMTrackerBase::FitTrack( const std::map<int,int> &hitcombo, double &xtrack, double &ytrack, double &xptrack, double &yptrack, double &chi2ndf, vector<double> &uresid, vector<double> &vresid ){

  //calculation of the best-fit line through the points was moved to its own function, since sometimes we want to perform ONLY that step; e.g.,
  //when calculating "exclusive" residuals:
  CalcLineOfBestFit( hitcombo, xtrack, ytrack, xptrack, yptrack );

  double chi2 = 0.0; 

  uresid.clear();
  vresid.clear();
  
  //I see no particularly good way to avoid looping over the hits again for the chi2 calculation:

  for( auto ilayer=hitcombo.begin(); ilayer != hitcombo.end(); ++ilayer ){
    int layer = ilayer->first;
    int hitidx = ilayer->second;

    int module = modindexhit2D[layer][hitidx];
    int clustidx = clustindexhit2D[layer][hitidx];

    double uhit = fModules[module]->fHits[clustidx].uhit;
    double vhit = fModules[module]->fHits[clustidx].vhit;

    TVector3 TrackOrigin( xtrack, ytrack, 0.0 );
    TVector3 TrackDirection( xptrack, yptrack, 1.0 );
    TrackDirection = TrackDirection.Unit();

    TVector2 UVtrack = GetUVTrack( module, TrackOrigin, TrackDirection ); 

    double utrack = UVtrack.X();
    double vtrack = UVtrack.Y();

    uresid.push_back( uhit - utrack );
    vresid.push_back( vhit - vtrack );
    
    chi2 += pow( (uhit-utrack)/fSigma_hitpos, 2 ) + pow( (vhit-vtrack)/fSigma_hitpos, 2 );
  }

  double ndf = double(2*nhits - 4);

  chi2ndf = chi2/ndf;
  
}

// "Odometer" algorithm for looping over possible combinations of one hit per layer:
bool SBSGEMTrackerBase::GetNextCombo( const std::set<int> &layers, const std::map<int,std::vector<int> > &hitlist, std::map<int,int> &hitcounter, std::map<int,int> &hitcombo, bool &firstcombo ){
  std::set<int>::iterator nextlayercounter = layers.begin();

  bool comboexists = true;
  
  for( auto layercounter = layers.begin(); layercounter != layers.end(); ++layercounter ){
    int layer = *layercounter;
    int nextlayer = *nextlayercounter;
    if( layer == nextlayer && !firstcombo ){
      if( hitcounter[layer]+1 < hitlist[layer].size() ){
	hitcounter[layer]++;
      } else {
	//reached last hit in current layer; roll back to first hit in this layer and increment hit counter in next layer:
	// Note that this means on the next iteration of layercounter,
	// layer == nextlayer will evaluate to true, and we will attempt to increment freehitcounter
	// for that layer as long as another free hit is available. Meanwhile, since the
	// free hit counter in the current layer has rolled back to 0, the next time we try to
	// populate a unique hit combination, starting from the first layer, we will loop over the hits
	// in the first layer again, having incremented the free hit counter for the next layer.
	// This process will repeat itself until we reach the last hit
	// in the last layer, at which point nextlayercounter will evaluate to layers.end,
	// and the iteration will stop.
	hitcounter[layer] = 0;
	++nextlayercounter;
      }
    }

    hitcombo[layer] = hitlist[layer][hitcounter[layer]];
    
    if( nextlayercounter == layers.end() ) comboexists = false; //we reached the last hit in the last layer. stop iteration
    
  }

  if( firstcombo ) firstcombo = false;
  
  return comboexists;
}

TVector3 SBSGEMTrackerBase::GetHitPosGlobal( int module, int clustindex ){
  
  //check that module and cluster are in range: these conditions should never evaluate to true, but we want to prevent
  // seg. faults anyway
  if( module < 0 || module >= fModules.size() ) return TVector3(-1.e12, -1.e12, -1.e12);
  if( clustindex < 0 || clustindex >= fModules[module]->fHits.size() ) return TVector3(-1.e12, -1.e12, -1.e12);
  
  return TVector3( fModules[module]->fHits[clustindex].xghit,
		   fModules[module]->fHits[clustindex].yghit,
		   fModules[module]->fHits[clustindex].zghit );
}

void SBSGEMTrackerBase::AddTrack( const std::map<int,int> &hitcombo, const vector<double> &BestTrack, double chi2ndf, const std::vector<double> &uresidbest, const std::vector<double> &vresidbest ){
  //AddTrack stores the best track found on each track-finding iteration in the appropriate data members of the class. It also takes care of
  //marking the hits on the track as used, and also marking all the 2D hits as used that contain any of the same 1D clusters as the found track:
  
  fNhitsOnTrack.push_back( hitcombo.size() );

  fXtrack.push_back( BestTrack[0] );
  fYtrack.push_back( BestTrack[1] );
  fXptrack.push_back( BestTrack[2] );
  fYptrack.push_back( BestTrack[3] );
  fChi2Track.push_back( chi2ndf );

  fresidu_hits.push_back( uresidbest );
  fresidv_hits.push_back( vresidbest );

  std::vector<int> modlisttemp,hitlisttemp;
  //We need to figure out a minimally expensive way to calculate "exclusive" residuals:

  //temporary vectors to hold exclusive residuals (residuals of the hit in question with respect to the track fitted to all the OTHER hits, excluding the hit in question):
  vector<double> eresidu, eresidv;
  
  for ( auto ilayer = hitcombo.begin(); ilayer != hitcombo.end(); ++ilayer ){
    int layer = ilayer->first;
    int hitidx = ilayer->second;

    int module = modindexhit2D[layer][hitidx];
    int iclust = clustindexhit2D[layer][hitidx];

    //Also: this is the time to mark the hits as used:
    hitused2D[layer][hitidx] = true;
    fModules[module]->fHits[iclust].ontrack = true;
    fModules[module]->fHits[iclust].trackidx = fNtracks_found;
    
    modlisttemp.push_back( module );
    hitlisttemp.push_back( iclust );

    //For exclusive residual calculation, we use CalcLineOfBestFit instead of FitTrack with the hit combo excluding the current layer:
    //copy the hit combo to a temporary local container:
    std::map<int,int> hitcombotemp = hitcombo;
    hitcombotemp.erase( layer ); //remove the current layer from the temporary copy of the list of hits

    //dummy variables to hold temporary track parameters:
    double xtemp,ytemp, xptemp,yptemp;

    // Calculate the line of best fit to all hits OTHER than the current one (without calculating chi2 or individual hit residuals):
    CalcLineOfBestFit( hitcombotemp, xtemp, ytemp, xptemp, yptemp );

    TVector3 TrackOrigin( xtemp, ytemp, 0.0 );
    TVector3 TrackDirection( xptemp, yptemp, 1.0 );
    TrackDirection = TrackDirection.Unit();

    TVector2 UVtrack = GetUVTrack( module, TrackOrigin, TrackDirection );

    TVector2 UVhit( fModules[module]->fHits[iclust].uhit, fModules[module]->fHits[iclust].vhit );

    eresidu.push_back( UVhit.X() - UVtrack.X() );
    eresidv.push_back( UVhit.Y() - UVtrack.Y() );
		      
  }
  
  fModListTrack.push_back( modlisttemp );
  fHitListTrack.push_back( hitlisttemp );

  feresidu_hits.push_back( eresidu );
  feresidv_hits.push_back( eresidv );

  //Purge hits containing either the same X cluster or the same Y cluster as any of the hits added to the track:
  //Note: This must be called AFTER adding the list of modules and the list of hits on the track to the track arrays:
  PurgeHits(fNtracks_found);
  
  fNtracks_found++;
}

//The purpose of this routine is to mark all the 2D hits as used that contain any of the same 1D U/V clusters as the 2D hits on this track.
//This routine accesses both the module hit arrays and the "hit list" arrays used by the track-finding. The loop over the entire 2D hit list in each layer
//on the track at the end of each track-finding iteration is maybe a bit expensive, but still probably cheap compared to the alternative (and will lead to fewer false tracks)
// The routine is designed to prevent re-use of the same 1D cluster in multiple tracks:
void SBSGEMTrackerBase::PurgeHits( int itrack ){
  for( int ihit=0; ihit<fNhitsOnTrack[itrack]; ihit++ ){
    int module = fModListTrack[itrack][ihit];
    int cluster = fHitListTrack[itrack][ihit];
    int layer = fModules[module]->fLayer;
    
    int uidx = fModules[module]->fHits[cluster].iuclust;
    int vidx = fModules[module]->fHits[cluster].ivclust;

    // Next we need to loop on the 2D hit list of this layer (the one used for track-finding) and mark any unused hits containing the same 1D (U/V) clusters
    // as the hits on this track as used:
    for( int jhit = 0; jhit<N2Dhits_layer[layer]; jhit++ ){
      int modj = modindexhit2D[layer][jhit];
      int clustj = clustindexhit2D[layer][jhit];

      if( modj == module ){
	int uj = fModules[modj]->fHits[clustj].iuclust;
	int vj = fModules[modj]->fHits[clustj].ivclust;

	if( uj == uidx || vj == vidx ){
	  hitused2D[layer][jhit] = true;
	}
      }
    }
  }
}

TVector3 SBSGEMTrackerBase::TrackIntersect( int module, TVector3 track_origin, TVector3 track_direction, double &sintersect ){
  TVector3 modpos = fModules[module]->GetOrigin();
  TVector3 modzaxis = fModules[module]->GetZax();

  sintersect = modzaxis.Dot( modpos - track_origin )/modzaxis.Dot( track_direction );

  return track_origin + sintersect * track_direction;
}

TVector2 SBSGEMTrackerBase::GetUVTrack( int module, TVector3 track_origin, TVector3 track_direction ){

  double sdummy; //we have to pass a double argument to hold the distance from origin to intersection:
  TVector3 TrackIntersect_Global = TrackIntersect( module, track_origin, track_direction, sdummy );
  TVector3 TrackIntersect_Local = fModules[module]->TrackToDetCoord( TrackIntersect_Global );

  TVector2 XYtrack( TrackIntersect_Local.X(), TrackIntersect_Local.Y() );
  return fModules[module]->XYtoUV( XYtrack );
}
