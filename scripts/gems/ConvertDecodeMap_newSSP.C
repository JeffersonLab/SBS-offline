#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <map>
#include <vector>
#include <set>

using namespace std;

void ConvertDecodeMap_newSSP( const char *oldfilename="gem_map_backup_July2021.txt", const char *newfilename="decode_map_for_analyzer_test.dat",
			      const char *detname="bb.gem"){

  ifstream infile(oldfilename);
  ofstream outfile(newfilename);

  TString currentline;


  //It appears that the only modification we need for the SSP version (Holly/Zeke/Sean's) is that the Crate ID is extracted from a dedicated line: 
  

  //let's make vectors mapped by layer and then by gempos:
  map<int,map<int,vector<int> > > crate, slot, prodID, axis, adcid, i2c, pos, invert, backplane, fiber;

  int lastcrate=-1;

  
  while( currentline.ReadLine(infile) ){
    if( currentline.Contains("Crate_") && currentline.Contains("Layer_") ){
      //extract Crate ID from 

      int istart = currentline.Index("Crate_") + 6;
      int istop = currentline.Index("Layer_");

      int len = istop-istart+1;

      TString sCrateID( currentline(istart,len) );

      sCrateID.ReplaceAll(" ","");
      //temporary hack to get things working; add 1 to lastcrate:
      lastcrate = sCrateID.Atoi() + 1;

      lastcrate = 20;
      std::cout << "found new crate ID line, Crate ID = " << lastcrate << endl;
    }
    
    if( currentline.BeginsWith("APV") ){ //every map entry starts with "APV";

      TObjArray *tokens = currentline.Tokenize(",");
      
      if( tokens->GetEntries() >= 13 ){ //valid entry

	TString slayer = ( (TObjString*) (*tokens)[2] )->GetString();
	int layer = slayer.Atoi();

	TString sgempos = ( (TObjString*) (*tokens)[12] )->GetString();
	int gempos = sgempos.Atoi();
	
	//TString scrate = ( (TObjString*) (*tokens)[1] )->GetString();
	

	TString sslot = ( (TObjString*) (*tokens)[1] )->GetString();

	TString sfiber = ( (TObjString*) (*tokens)[3] )->GetString();

	TString sprodID = ( (TObjString*) (*tokens)[4] )->GetString();
	

	TString saxis = ( (TObjString*) (*tokens)[5] )->GetString();
	

	TString sadcid = ( (TObjString*) (*tokens)[6] )->GetString();
	

	TString si2c = ( (TObjString*) (*tokens)[7] )->GetString();
	

	TString spos = ( (TObjString*) (*tokens)[8] )->GetString();

	TString sinvert = ( (TObjString*) (*tokens)[9] )->GetString();
	

	if( lastcrate >= 0 ){
	
	  crate[layer][gempos].push_back( lastcrate );
	  //slot[layer][gempos].push_back( sslot.Atoi() );
	  slot[layer][gempos].push_back( 11 );
	  fiber[layer][gempos].push_back( sfiber.Atoi() );
	  prodID[layer][gempos].push_back( sprodID.Atoi() );
	  axis[layer][gempos].push_back( saxis.Atoi() );
	  adcid[layer][gempos].push_back( sadcid.Atoi() );
	  i2c[layer][gempos].push_back( si2c.Atoi() );
	  pos[layer][gempos].push_back( spos.Atoi() );
	  invert[layer][gempos].push_back( sinvert.Atoi() );
	}
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
      modentry.Form( "%s.m%d.chanmap = ", detname, modcounter );

      outfile << modentry << endl;
      
      for( int mapentry=0; mapentry<gemmap[module].size(); mapentry++ ){
	TString smapline;

	//Slot is used twice in a redundant fashion. We probably want to change that later to avoid confusion:
	smapline.Form( " %6d  %6d  %6d  %6d  %6d  %6d  %6d  %6d  %6d ",
		       crate[layer][module][mapentry],
		       slot[layer][module][mapentry],
		       fiber[layer][module][mapentry],
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
