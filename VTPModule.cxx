/*
 SBS HCAL/ECAL VTP data format
 Based on the manual from B. Raydo (Jan. 28, 2025)
 Modified from NPSlib/VTPModule.cxx 

  Data Type list  
  0 Block header
  1 Block Trailer
  2 Event Header
  3 Trigger Time
  12 Expanded (Data SubType)
  12.11 HCAL Cluster
  12.12 ECAL Cluster
  
  Data words: type (bit 31) + 4 bits for type tag (30-27)
  Type defining word 31 = 1, continuation 31 = 0
  
  ------------------
  Block Header (type 0)
  31    = 1
  30-27 = 0
  26-22 = SLOTID
  21-18 = UNDEFINED
  17-8  = EVENT_PER_BLOCK
  7 -0  = BLOCK_CNT (Event block number)
  
  Block Trailer (type 1)
  31    = 1
  30-27 = 1
  26-22 = SLOTID
  21-0  = NUM_WORDS
  17-8  = EVENT_PER_BLOCK
  7 -0  = BLOCK_CNT (Event block number)
  
  Event Header(type 2)
  31    = 1
  30-27 = 2
  26-22 = UNDEFINED
  21-0  = TRIGGER_NUMBER

  Trigger Time (type 3)  - 2 words
  Word 1:
  31    = 1
  30-27 = 3
  26-24 = UNDEFINED
  23-0  = TRIGGER_TIME_L (lower 24 bits)

  Word 2:  
  31    = 0
  30-24 = UNDEFINED
  23-0  = TRIGGER_TIME_H (upper 24 bits)

  SBS HCAL Cluster (type 12) - 2 words
  Word 1:
  31    = 1
  30-27 = 12
  26-23 = 11
  22-16 = -
  15-0  = E (16bits cluster energy)

  Word 2:
  31    = 0
  30-24 = -
  23-20 = 4 bits cluster column coordinate 
  19-15 = 5 bits cluster row coordinate
  14-11 = 4 bits number of hits in the cluster
  10-0  = 11 bits cluster time in 4ns uinits

  SBS ECAL Cluster (type 12) - 2 words
  Word 1:
  31    = 1
  30-27 = 12
  26-23 = 12
  22-14 = -
  13-0  = E (14bits cluster energy)

  Word 2:
  31    = 0
  30-27 = -
  26-22 = 5 bits cluster column coordinate 
  21-15 = 7 bits cluster row coordinate
  14-11 = 4 bits number of hits in the cluster
  10-0  = 11 bits cluster time in 4ns uinits
*/

#include "VTPModule.h"
#include "THaSlotData.h"
#include "TMath.h"

#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <cassert>
#include <stdexcept>
#include <map>
#include <sstream>

using namespace std;

#ifdef DEBUG
#include <fstream>
#endif

namespace Decoder {

using vsiz_t = vector<int>::size_type;

Module::TypeIter_t VTPModule::fgThisType =
        DoRegister( ModuleType( "Decoder::VTPModule" , 514 ));

//_____________________________________________________________________________
VTPModule::VTPModule(UInt_t crate, UInt_t slot)
  : PipeliningModule(crate, slot), 
    vtp_header_data{}, vtp_trigger_data{}, vtp_cluster_data{},
    block_header_found(false), block_trailer_found(false),
    event_header_found(false), slots_match(false)
{
  IsInit = false;
  VTPModule::Init();
}

//_____________________________________________________________________________
VTPModule::~VTPModule() = default;

//_____________________________________________________________________________
void VTPModule::Clear( Option_t* opt )
{
  // Clear event-by-event data
  PipeliningModule::Clear(opt);
  vtp_header_data.clear();
  vtp_trigger_data.clear();
  vtp_cluster_data.clear();
  
  // Initialize data_type_def to FILLER and data types to false
  data_type_def = 15;
  // Initialize data types to false
  block_header_found = block_trailer_found = event_header_found = slots_match = false;
}

//_____________________________________________________________________________
void VTPModule::Init()
{
  Module::Init();
#if defined DEBUG && defined WITH_DEBUG
  // This will make a HUGE output
  delete fDebugFile; fDebugFile = 0;
  fDebugFile = new ofstream;
  fDebugFile->open(string("VTP_debug.txt"));
#endif
  //fDebugFile=0;
  Clear();
  IsInit = kTRUE;
  fName = "VTP";
}

//_____________________________________________________________________________
  Int_t VTPModule::Decode( const UInt_t* pdat )
{
  assert(pdat);
  uint32_t data = *pdat;
  uint32_t data_type_id = (data >> 31) & 0x1;  // Data type identification, mask 1 bit
  if( data_type_id == 1 )
    data_type_def = (data >> 27) & 0xF;        // Data type defining words, mask 4 bits

  // Debug output
#ifdef WITH_DEBUG
  if( fDebugFile )
    *fDebugFile << "VTPModule::Decode:: VTP DATA TYPES"
                << " >> data = " << hex << data << dec
                << " >> data word id = " << data_type_id
                << " >> data type = " << data_type_def
                << endl;

#endif

  // Ensure that slots match and do not decode if they differ.
  // This should never happen if PipeliningModule::LoadBank selected the
  // correct bank.
  if( !slots_match && data_type_def != 0 ) {
#ifdef WITH_DEBUG
    if( fDebugFile )
      *fDebugFile << "VTPModule::Decode:: "
                  << "fSlot & VTP slot do not match AND data type != 0"
                  << endl;
#endif
    return -1;
  }
  
  // Acquire data objects depending on the data type defining word
  switch( data_type_def ) {
    case 0: // Block header, indicates the beginning of a block of events
      DecodeBlockHeader(data, data_type_id);
      break;
    case 1: // Block trailer, indicates the end of a block of events
      DecodeBlockTrailer(data);
      break;
    case 2: // Event header, indicates start of an event, includes the trigger number
      DecodeEventHeader(data);
      break;
    case 3:  // Trigger time, time of trigger occurrence relative to the most recent global reset
      DecodeTriggerTime(data, data_type_id);
      break;
      break;
    case 12: // cluster
      DecodeCluster(data, data_type_id);
      break;
    case 13: // trigger decision
      DecodeTriggerDecision(data, data_type_id);
      break;
    case 4:  // Undefined type
    case 5:  // Undefined type
    case 6:  // Undefined type
    case 7:  // Undefined type
    case 8:  // Undefined type
    case 9:  // Undefined type
    case 10: // Undefined type
    case 11: // Undefined type
    case 14: // Data not valid
    case 15: // Filler Word, should be ignored
      UnsupportedType(data, data_type_id);
      break;
  default:
    throw logic_error("VTPModule: incorrect masking of data_type_def");
  }  // data_type_def switch
  
#ifdef WITH_DEBUG
  if( fDebugFile )
    *fDebugFile << "**********************************************************************"
                << "\n" << endl;
#endif
  return block_trailer_found;
}

//_____________________________________________________________________________
void VTPModule::DecodeBlockHeader( UInt_t pdat, uint32_t data_type_id )
{
  static const string here{ "VTPModule::DecodeBlockHeader" };

  if( data_type_id ) {
    block_header_found = true;                     // Set to true if found
    vtp_header_data.slot_blk_hdr = (pdat >> 22) & 0x1F;  // Slot number (set by VME64x backplane), mask 5 bits
    // Debug output
#ifdef WITH_DEBUG
    if( fDebugFile )
      *fDebugFile << "VTPModule::Decode:: Slot from VTP block header = " << vtp_header_data.slot_blk_hdr << endl;
#endif
    // Ensure that slots from cratemap and VTP match
    slots_match = ((uint32_t) fSlot == vtp_header_data.slot_blk_hdr);
    if( !slots_match )
      return;
    // vtp_header_data.iblock_num = (pdat >> 8) & 0x3FF;    // Event block number, mask 10 bits
    // vtp_header_data.nblock_events = (pdat >> 0) & 0xFF;  // Number of events in block, mask 8 bits
    vtp_header_data.nblock_events = (pdat >> 0) & 0x3FF; // Number of events in block, mask 10 bits
    vtp_header_data.iblock_num = (pdat >> 8) & 0xFF;    // Event block number, mask 8 bits

    // Debug output
#ifdef WITH_DEBUG
    if( fDebugFile )
      *fDebugFile << "VTPModule::Decode:: VTP BLOCK HEADER"
                  << " >> data = " << hex << pdat << dec
                  << " >> slot = " << vtp_header_data.slot_blk_hdr
                  << " >> event block number = " << vtp_header_data.iblock_num
                  << " >> num events in block = " << vtp_header_data.nblock_events
                  << endl;
#endif
  }
}

//_____________________________________________________________________________
void VTPModule::DecodeBlockTrailer( UInt_t pdat )
{
  block_trailer_found = true;
  vtp_header_data.slot_blk_trl = (pdat >> 22) & 0x1F;       // Slot number (set by VME64x backplane), mask 5 bits
  vtp_header_data.nwords_inblock = (pdat >> 0) & 0x3FFFFF;  // Total number of words in block of events, mask 22 bits
  // Debug output
#ifdef WITH_DEBUG
  if( fDebugFile )
    *fDebugFile << "VTPModule::Decode:: VTP BLOCK TRAILER"
                << " >> data = " << hex << pdat << dec
                << " >> slot = " << vtp_header_data.slot_blk_trl
                << " >> nwords in block = " << vtp_header_data.nwords_inblock
                << endl;
#endif
}

//_____________________________________________________________________________
void VTPModule::DecodeEventHeader( UInt_t pdat )
{
  event_header_found = true;
  vtp_header_data.trig_num = (pdat >> 0) & 0x3FFFFF;  // Trigger number, mask 22 bits
#ifdef WITH_DEBUG
  if( fDebugFile )
    *fDebugFile << "VTPModule::Decode:: VTP EVENT HEADER"
                << " >> data = " << hex << pdat << dec
                << " >> trigger number = " << vtp_header_data.trig_num
                << endl;
#endif
}

//_____________________________________________________________________________
void VTPModule::DecodeTriggerTime( UInt_t pdat, uint32_t data_type_id )
{
  if( data_type_id )  // Trigger time word 1
    vtp_header_data.trig_time_w1 = (pdat >> 0) & 0xFFFFFF;  // Time = T_D T_E T_F (lower 24bits)
  else  // Trigger time word 2
    vtp_header_data.trig_time_w2 = (pdat >> 0) & 0xFFFFFF;  // Time = T_A T_B T_C (upper 24bits)
  // Time = T_A T_B T_C T_D T_E T_F
  vtp_header_data.trig_time = (vtp_header_data.trig_time_w2 << 24) | vtp_header_data.trig_time_w1;
  // Debug output
#ifdef WITH_DEBUG
  if( fDebugFile )
    *fDebugFile << "VTPModule::Decode:: VTP TRIGGER TIME"
                << " >> data = " << hex << pdat << dec
                << " >> trigger time word 1 = " << vtp_header_data.trig_time_w1
                << " >> trigger time word 2 = " << vtp_header_data.trig_time_w2
                << " >> trigger time = " << vtp_header_data.trig_time
                << endl;
#endif
}
  
//_____________________________________________________________________________
void VTPModule::DecodeCluster( UInt_t pdat, uint32_t data_type_id )
{
  uint32_t ce = 0;
  uint32_t ct = 0;
  uint32_t cn = 0;
  uint32_t cx = 0;
  uint32_t cy = 0;
  if( data_type_id )  { //  cluster word 1
    uint32_t det_id = (pdat >> 23) & 0xF; // 12: ECAL 13: HCAL
    vtp_cluster_data.det_id = det_id;
      ce = (pdat >> 0) & 0x3FFF; //  cluster energy, mask 14 bits
    vtp_cluster_data.energy.push_back( ce ); 
    
#ifdef WITH_DEBUG
    if( fDebugFile )
      *fDebugFile << "VTPModule::Decode:: VTP  Cluster"
                  << " >> data = " << hex << pdat << dec
                  << " >> energy = " << ce << endl;
#endif

  }
  else { //  ckuster word 2
 
    if(vtp_cluster_data.det_id == HCAL) {
      cx = (pdat >> 20) & 0xF;    //  cluster x (col) coordinate, mask 4 bits
      cy = (pdat >> 15) & 0x1F;   //  cluster y (row) coordinate, mask 5 bits
      cn = (pdat >> 11) & 0xF;    //  number of hits in the cluster, mask 4 bits
      ct = (pdat >> 0)  & 0x7FF;  //  cluster time, mask 11 bits
    }
    else {
      cx = (pdat >> 22) & 0x1F;   //  cluster x coordinate, mask 5 bits
      cy = (pdat >> 15) & 0x7F;   //  cluster y coordinate, mask 7 bits
      cn = (pdat >> 11) & 0xF;    //  cluster n blocks, mask 4 bits
      ct = (pdat >> 0)  & 0x7FF;  //  cluster time, mask 11 bits
    }

    vtp_cluster_data.time.push_back( ct );
    vtp_cluster_data.nblocks.push_back( cn );
    vtp_cluster_data.xcoord.push_back( cx );
    vtp_cluster_data.ycoord.push_back( cy );

#ifdef WITH_DEBUG
    if( fDebugFile )
      *fDebugFile << "VTPModule::Decode:: VTP  Cluster"
                  << " >> n blocks = " << cn
                  << " >> time  = " << ct
                  << " >> x coordinate = " << cx
                  << " >> y coordinate = " << cy
                  << endl;
#endif
  } 
}
  
//_____________________________________________________________________________
  void VTPModule::DecodeTriggerDecision( UInt_t pdat, uint32_t data_type_id )
  {
    uint32_t tpattern = 0;
    uint32_t ttime    = 0;
    uint32_t ttype0   = 0;
    uint32_t ttype1   = 0;
    uint32_t ttype2   = 0;
    uint32_t ttype3   = 0;
    uint32_t ttype4   = 0;
    uint32_t ttype5   = 0;
    
    if( data_type_id ) {  // trigger decision word 1
      
      tpattern = (pdat >> 0) & 0xFFFF; // trigger pattern, mask 16 bits
      ttime    = (pdat >> 16) & 0x7F;  // trigger time, mask 11 bits
      
      if( (tpattern >> 0 ) & 0x1 ) ttype0 = 1;
      if( (tpattern >> 1 ) & 0x1 ) ttype1 = 1;
      if( (tpattern >> 2 ) & 0x1 ) ttype2 = 1;
      if( (tpattern >> 3 ) & 0x1 ) ttype3 = 1;
      if( (tpattern >> 4 ) & 0x1 ) ttype4 = 1;
      if( (tpattern >> 5 ) & 0x1 ) ttype5 = 1;

      //      cout << " tpattern = " << tpattern << endl;
      vtp_trigger_data.trigtype0.push_back( ttype0 );
      vtp_trigger_data.trigtype1.push_back( ttype1 );
      vtp_trigger_data.trigtype2.push_back( ttype2 );
      vtp_trigger_data.trigtype3.push_back( ttype3 );
      vtp_trigger_data.trigtype4.push_back( ttype4 );
      vtp_trigger_data.trigtype5.push_back( ttype5 );
      vtp_trigger_data.trigtime.push_back( ttime );
    }
    
#ifdef WITH_DEBUG
    if( fDebugFile )
      *fDebugFile << "VTPModule::Decode:: VTP Trigger"
		  << " >> data  = " << hex << pdat << dec
		  << " >> trig0 = " << ttype0
		  << " >> trig1 = " << ttype1
		  << " >> trig2 = " << ttype2
		  << " >> trig time = " << ttime
		  << endl;
#endif
  }

//_____________________________________________________________________________
void VTPModule::UnsupportedType( UInt_t pdat, uint32_t data_type_id )
{
  // Handle unsupported, invalid, or irrelevant/non-decodable data types
#ifdef WITH_DEBUG
  // Data type descriptions
  static const vector<string> what_text{ "UNDEFINED TYPE",
                                         "DATA NOT VALID",
                                         "FILLER WORD",
                                         "INCORRECT DECODING" };
  // Lookup table data_type -> message number
  static const map<uint32_t, uint32_t> what_map = {
    // undefined type
    { 4,  0 },
    { 5,  0 },
    { 6,  0 },
    { 7,  0 },
    { 8,  0 },
    { 9,  0 },
    { 10, 0 },
    { 11, 0 },
    // data not valid
    { 14, 1 },
    // filler word
    { 15, 2 }
  };
  auto idx = what_map.find(data_type_def);
  // Message index. The last message means this function was called when
  // it shouldn't have been called, i.e. coding error in DecodeOneWord
  size_t i = (idx == what_map.end()) ? what_text.size() - 1 : idx->second;
  const string& what = what_text[i];
  ostringstream str;
  str << "VTPModule::Decode:: " << what
      << " >> data = " << hex << pdat << dec
      << " >> data type id = " << data_type_id
      << endl;
  if( fDebugFile )
    *fDebugFile << str.str();
  if( idx == what_map.end() )
    cerr << str.str();
#endif
}

//_____________________________________________________________________________
UInt_t VTPModule::GetHeaderTriggerTime() 
{
  // Truncate to 32 bits
  UInt_t shorttime = vtp_header_data.trig_time;
#ifdef WITH_DEBUG
  if( fDebugFile )
    *fDebugFile << "VTPModule::GetTriggerTime = "
                << vtp_header_data.trig_time << " " << shorttime << endl;
#endif
  return shorttime;
}

//_____________________________________________________________________________
void VTPModule::LoadTHaSlotDataObj( THaSlotData* sldat )
{
  // Load THaSlotData

  // trigger data
  for( vsiz_t itrig = 0; itrig < vtp_trigger_data.trigtype0.size(); itrig++ )
    sldat->loadData("scaler", 0, vtp_trigger_data.trigtype0[itrig], vtp_trigger_data.trigtype0[itrig]);
  for( vsiz_t itrig = 0; itrig < vtp_trigger_data.trigtype1.size(); itrig++ )
    sldat->loadData("scaler", 0, vtp_trigger_data.trigtype1[itrig], vtp_trigger_data.trigtype1[itrig]);
  for( vsiz_t itrig = 0; itrig < vtp_trigger_data.trigtype2.size(); itrig++ )
    sldat->loadData("scaler", 0, vtp_trigger_data.trigtype2[itrig], vtp_trigger_data.trigtype2[itrig]);
  for( vsiz_t itrig = 0; itrig < vtp_trigger_data.trigtime.size(); itrig++ )
    sldat->loadData("scaler", 0, vtp_trigger_data.trigtime[itrig], vtp_trigger_data.trigtime[itrig]);

  // cluster data
  sldat->loadData("scaler", 0, vtp_cluster_data.det_id, vtp_cluster_data.det_id);
  for( vsiz_t iclus = 0; iclus < vtp_cluster_data.energy.size(); iclus++ )
    sldat->loadData("scaler", 0, vtp_cluster_data.energy[iclus], vtp_cluster_data.energy[iclus]);
  for( vsiz_t iclus = 0; iclus < vtp_cluster_data.time.size(); iclus++ )
    sldat->loadData("scaler", 0, vtp_cluster_data.time[iclus], vtp_cluster_data.time[iclus]);
  for( vsiz_t iclus = 0; iclus < vtp_cluster_data.nblocks.size(); iclus++ )
    sldat->loadData("scaler", 0, vtp_cluster_data.nblocks[iclus], vtp_cluster_data.nblocks[iclus]);
  for( vsiz_t iclus = 0; iclus < vtp_cluster_data.xcoord.size(); iclus++ )
    sldat->loadData("scaler", 0, vtp_cluster_data.xcoord[iclus], vtp_cluster_data.xcoord[iclus]);
  for( vsiz_t iclus = 0; iclus < vtp_cluster_data.ycoord.size(); iclus++ )
    sldat->loadData("scaler", 0, vtp_cluster_data.ycoord[iclus], vtp_cluster_data.ycoord[iclus]);

}

//_____________________________________________________________________________
UInt_t VTPModule::LoadSlot( THaSlotData* sldat, const UInt_t* evbuffer,
                                const UInt_t* pstop )
{
  // Load from evbuffer between [evbuffer,pstop]
  return LoadSlot(sldat, evbuffer, 0, pstop + 1 - evbuffer);
}

//_____________________________________________________________________________
UInt_t VTPModule::LoadSlot( THaSlotData *sldat, const UInt_t* evbuffer,
                                UInt_t pos, UInt_t len)
{
  // Load from bank data in evbuffer between [pos,pos+len)
  const auto* p = evbuffer + pos;
  const auto* q = p + len;
  while( p != q ) {
    if( Decode(p++) == 1 )
      break;  // block trailer found
  }

  LoadTHaSlotDataObj(sldat);
  
  return fWordsSeen = p - (evbuffer + pos);
}

//_____________________________________________________________________________


} // namespace Decoder

ClassImp(Decoder::VTPModule)
