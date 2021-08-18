#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <map>
#include <vector>
#include <set>

using namespace std;

void ConvertDecodeMap( const char *oldfilename="database/gem_map_uva_eel_3crates_layer12345.cfg", const char *newfilename="decode_map_for_analyzer.dat" ){

  ifstream infile(oldfilename);
  ofstream outfile(newfilename);

  TString currentline;

  //Siyu's version of the map file specifies layer and GEMPOS (position within a layer):
  

  //let's make vectors mapped by layer and then by gempos:
  map<int,map<int,vector<int> > > crate, slot, prodID, axis, adcid, i2c, pos, invert, backplane;
  
  while( currentline.ReadLine(infile) ){ 
    if( currentline.BeginsWith("APV") ){ //every map entry starts with "APV";

      TObjArray *tokens = currentline.Tokenize(",");
      
      if( tokens->GetEntries() >= 13 ){ //valid entry

	TString slayer = ( (TObjString*) (*tokens)[2] )->GetString();
	int layer = slayer.Atoi();

	TString sgempos = ( (TObjString*) (*tokens)[12] )->GetString();
	int gempos = sgempos.Atoi();
	
	TString scrate = ( (TObjString*) (*tokens)[1] )->GetString();
	

	TString sslot = ( (TObjString*) (*tokens)[3] )->GetString();


	TString sprodID = ( (TObjString*) (*tokens)[4] )->GetString();
	

	TString saxis = ( (TObjString*) (*tokens)[5] )->GetString();
	

	TString sadcid = ( (TObjString*) (*tokens)[6] )->GetString();
	

	TString si2c = ( (TObjString*) (*tokens)[7] )->GetString();
	

	TString spos = ( (TObjString*) (*tokens)[8] )->GetString();

	TString sinvert = ( (TObjString*) (*tokens)[9] )->GetString();
	

	crate[layer][gempos].push_back( scrate.Atoi() );
	slot[layer][gempos].push_back( sslot.Atoi() );
	prodID[layer][gempos].push_back( sprodID.Atoi() );
	axis[layer][gempos].push_back( saxis.Atoi() );
	adcid[layer][gempos].push_back( sadcid.Atoi() );
	i2c[layer][gempos].push_back( si2c.Atoi() );
	pos[layer][gempos].push_back( spos.Atoi() );
	invert[layer][gempos].push_back( sinvert.Atoi() );
      }      
    }
  }

  int modcounter = 0;
  
  for( auto ilayer = crate.begin(); ilayer != crate.end(); ++ilayer ){
    int layer = ilayer->first;
    map<int,vector<int> > gemmap = crate[layer];
    for( auto igem = gemmap.begin(); igem != gemmap.end(); ++igem ){
      int module = igem->first;
      // layer and module are sorted now:
      TString modentry;
      modentry.Form( "sbs.uvagem.m%d.chanmap = ", modcounter );

      outfile << modentry << endl;
      
      for( int mapentry=0; mapentry<gemmap[module].size(); mapentry++ ){
	TString smapline;

	//Slot is used twice in a redundant fashion. We probably want to change that later to avoid confusion:
	smapline.Form( " %6d  %6d  %6d  %6d  %6d  %6d  %6d  %6d  %6d ",
		       crate[layer][module][mapentry],
		       slot[layer][module][mapentry],
		       slot[layer][module][mapentry],
		       prodID[layer][module][mapentry],
		       adcid[layer][module][mapentry],
		       i2c[layer][module][mapentry],
		       pos[layer][module][mapentry],
		       invert[layer][module][mapentry],
		       axis[layer][module][mapentry] );

	outfile << smapline << endl;
      }
      

      modcounter++; 
    }
  }
  
}
