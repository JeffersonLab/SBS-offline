//*-- Author :   Paul King  13-April-2025
////////////////////////////////////////////////////////////////////////
//
// SBSHelicityDecoderBoard
//
////////////////////////////////////////////////////////////////////////

#include "SBSHelicityDecoderBoard.h"
#include "THaEvData.h"
#include "VarDef.h"
#include <iostream>

using namespace std;

const Int_t  SBSHelicityDecoderBoard::fNumDecoderWords = 14;  /// The number of data words in the HD bank; there are at least 4 header/trailer words also

//____________________________________________________________________
SBSHelicityDecoderBoard::SBSHelicityDecoderBoard( const char* name, const char* description,
			  THaApparatus* app ) : 
  THaHelicityDet(name, description, app),
  fHelicityDelay(0)
{

}

//____________________________________________________________________
SBSHelicityDecoderBoard::SBSHelicityDecoderBoard()
  : THaHelicityDet(),
    fHelicityDelay(0)
{
  // Default constructor - for ROOT I/O only
}

//____________________________________________________________________
SBSHelicityDecoderBoard::~SBSHelicityDecoderBoard() 
{
  // Destructor

  RemoveVariables();
}

//_____________________________________________________________________________
Int_t SBSHelicityDecoderBoard::Begin( THaRunBase* )
{
  return 0;
}

//_____________________________________________________________________________
Int_t SBSHelicityDecoderBoard::End( THaRunBase* )
{
   return 0;
}


//____________________________________________________________________
Int_t SBSHelicityDecoderBoard::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  // Initialize global variables

  // Define standard variables from base class
  Int_t ret = THaHelicityDet::DefineVariables( mode );
  if( ret )
    return ret;

  // Define variables
  const RVarDef var[] = {
    { "seed_reported",    "Seed for Reported Helicity",       "fSeed_Reported" },
    { "seed_true",        "Seed for True Helicity",           "fSeed_True" },
    { "n_tstable_rise",   "Number of Tstable rise signals",   "fNum_TStable_Fall" },
    { "n_tstable_fall",   "Number of Tstable fall signals",   "fNum_TStable_Rise" },
    { "n_patsync",        "Number of pattern sync signals",   "fNum_Pattern_Sync" },
    { "n_pairsync",       "Number of pair sync signals",      "fNum_Pair_Sync" },
    { "time1",            "Trigger time since Tstable start", "fTime_since_TStable" },
    { "time2",            "Trigger time since Tsettle start", "fTime_since_TSettle" },
    { "duration_tstable", "Trigger time since Tsettle start", "fLast_Duration_TStable" },
    { "duration_tsettle", "Trigger time since Tsettle start", "fLast_Duration_TSettle" },
    { nullptr }
  };
  return DefineVarsFromList( var, mode );
}

//____________________________________________________________________
Int_t SBSHelicityDecoderBoard::ReadDatabase( const TDatime& date )
{
  // Read database

  // Read ADCHelicity parameters
  Int_t st = THaHelicityDet::ReadDatabase( date );
  if( st != kOK )
    return st;

  // Read parameters for this class
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  DBRequest req[] = {
    { "verbose",       &fDebug,         kInt, 0, true, 1 },
    { "heldelay",      &fHelicityDelay, kInt, 0, true, 1},
    { "rocid",         &ROC_ID,         kInt, 0, false, 0},
    { "bankid",        &Bank_ID,        kInt, 0, false, 0},
    { nullptr }
  };
  st = LoadDB( file, date, req );
  fclose(file);
  if( st )
    return kInitError;

  return kOK;
}

//____________________________________________________________________
void SBSHelicityDecoderBoard::Clear( Option_t* opt ) 
{
  // Clear the event data
  THaHelicityDet::Clear(opt);

  //  HelicityDecoder data values
  fSeed_Reported               = 0;
  fNum_TStable_Fall            = 0;
  fNum_TStable_Rise            = 0;
  fNum_Pattern_Sync            = 0;
  fNum_Pair_Sync               = 0;
  fTime_since_TStable          = 0;
  fTime_since_TSettle          = 0;
  fLast_Duration_TStable       = 0;
  fLast_Duration_TSettle       = 0;
  fPatternPhase                = 0;
  fEventPolarity               = 0;
  fReportedPatternHel          = 0;
  fBit_Helicity                = 0;
  fBit_PairSync                = 0;
  fBit_PatSync                 = 0;
  fBit_TStable                 = 0;
  fEvtHistory_PatSync          = 0;
  fEvtHistory_PairSync         = 0;
  fEvtHistory_ReportedHelicity = 0;
  fPatHistory_ReportedHelicity = 0;

  
}




//____________________________________________________________________
Int_t SBSHelicityDecoderBoard::CheckPredictor(UInt_t& testval)
{
  Int_t  retval;
  UInt_t rval = (testval)&0x3fffffff;
  UInt_t lval = (testval>>2)&0x3fffffff;
  GetRandbit30(lval);
  GetRandbit30(lval);
  if (fDebug>0) std::cerr << std::hex << "Original seed: " << testval
			  << "; rval:  " << rval
			  << "; lval:  " << ((testval>>2)&0x3fffffff)
			  << "; regnerated seed: " << lval
			  << std::endl;
  if (rval == lval) {
    if (fDebug>0) std::cout << "Helicity seed is good" << std::endl;
    retval = kOK;
  } else {
    if (fDebug>0) std::cerr << "PROBLEM reconstructing helicity seed!!" << std::endl;
    retval = 10;
  }
  return retval;
}

//____________________________________________________________________
UInt_t SBSHelicityDecoderBoard::GetRandbit30(UInt_t& ranseed)
{
  /** This is a 30 bit random bit generator according to the "new" algorithm
      described in "G0 Helicity Digital Controls" by E. Stangland, R. Flood, H. Dong.


      The helicity board uses a maximum-length linear feedback shift registers
      for the generation of a pseudo-random sequence of bits.  The length of the
      register (24 bits or 30 bits) defines the length before a sequence is
      repeated: 2^n - 1.

      For a mathematical introduction to the generation of pseudo-random numbers
      with maximum-length linear feedback shift registers (LFSR), see the
      following web references:
         http://en.wikipedia.org/wiki/Linear_feedback_shift_register
	 http://www.newwaveinstruments.com/resources/articles/m_sequence_linear_feedback_shift_register_lfsr.htm

      In particular, the used solutions are to place XNOR taps at the bits
         24 stages, 4 taps:  (47 sets)
           [24, 23, 21, 20]
         30 stages, 4 taps:  (104 sets)
           [30, 29, 28, 7]
  **/

  UInt_t bit7    = (ranseed & 0x00000040) != 0;
  UInt_t bit28   = (ranseed & 0x08000000) != 0;
  UInt_t bit29   = (ranseed & 0x10000000) != 0;
  UInt_t bit30   = (ranseed & 0x20000000) != 0;

  UInt_t result = (bit30 ^ bit29 ^ bit28 ^ bit7) & 0x1;

  ranseed =  ( (ranseed << 1) | result ) & 0x3FFFFFFF;

  return(result);
}


//____________________________________________________________________
Int_t SBSHelicityDecoderBoard::Decode( const THaEvData& evdata )
{

  Int_t st = THaHelicityDet::Decode( evdata );
  if( st < 0 )
    return st;

  //  Read the Helicity Decoder data bank
  BANKinfo bankinfo;
  //  Set the bank info for the ring buffer entries
  bankinfo.roc = ROC_ID;
  bankinfo.bankid = Bank_ID;
  const UInt_t *lbuff;
  lbuff = GetBankBuffer(evdata, bankinfo);
  if ( ! DecodeBank(lbuff, bankinfo)){
    return (-1);
  } else {
    //  We have correctly filled the data members for this event.
    //  Use the history of the last 32 patterns to cross-check
    //  the predictor, then calulate the true helicity.
    if (CheckPredictor(fPatHistory_ReportedHelicity)==kOK){
      fSeed_True=fSeed_Reported&0x3fffffff;
      if (fHelicityDelay>0){
	for (size_t ipat=0; ipat<fHelicityDelay; ipat++){
	  GetRandbit30(fSeed_True);
	}
      }
      Int_t event_true_hel = (fSeed_True&1) ^ fEventPolarity;
      if (fBit_TStable==1){
	if ((event_true_hel)==1){
	  fHelicity = kPlus;
	} else {
	  fHelicity = kMinus;
	}
      } else {
	//  Mark the helicity as unknown during the spin-flip
	fHelicity = kUnknown;
      }
    } else {
      fSeed_True = 0;
      fHelicity = kUnknown;
      std::cerr << "HELICTY PREDICTION FAILURE" << std::endl;
    }
  }
  return kOK;
}

//____________________________________________________________________
const UInt_t* SBSHelicityDecoderBoard::GetBankBuffer( const THaEvData& evdata,
      BANKinfo& info )
{
   UInt_t loff; // Local offset
   UInt_t banklen, bankhead;
   UInt_t len = evdata.GetRocLength(info.roc);
   if (len <= 4) {
      info.length = 0;
      return nullptr;
   }
   const UInt_t* buff = evdata.GetRawDataBuffer(info.roc);
   UInt_t rochead = buff[1];
   if ((rochead&0xff00)!=0x1000){
      std::cerr <<  "SBSHelicityDecoderBoard::GetBankBuffer:  ROC "
		<<  std::dec << info.roc
		<<  " is not made of banks"  <<std::endl;
      info.length = 0;
      return nullptr;
   }
   loff = 2;
   while (loff < len-2){
      banklen  = buff[loff];
      bankhead = buff[loff+1];
      if (((bankhead&0xffff0000)>>16)==info.bankid) break;
      loff += 1+banklen;
   }
   if ( (((bankhead&0xffff0000)>>16)==info.bankid) &&
         ((loff+1+banklen) <= len+1) ){
      buff+=loff+2;
      info.length = banklen-1;
   } else {
      info.length = 0;
      return nullptr;
   }
   return buff;
}

//____________________________________________________________________
Bool_t  SBSHelicityDecoderBoard::DecodeBank(const UInt_t *lbuff, const BANKinfo& info)
{
  uint32_t type_last = 15;       /* initialize to type FILLER WORD */
  uint32_t time_last = 0;
  uint32_t decoder_index = 0;
  uint32_t num_decoder_words = 1;
  
  uint32_t slot_id_ev_hd = 0;
  uint32_t slot_id_dnv = 0;
  uint32_t slot_id_fill = 0;
   
  uint32_t new_type;
  uint32_t type;
  uint32_t slot_id_hd;
  uint32_t mod_id_hd;
  uint32_t slot_id_tr;
  uint32_t n_evts;
  uint32_t blk_num;
  uint32_t n_words;
  uint32_t evt_num_1;
  uint32_t trig_time;
  uint32_t time_now;
  uint32_t time_1;
  uint32_t time_2;
  uint32_t num_words;

  uint32_t data;

  if (info.length<fNumDecoderWords+4){
    std::cerr << "Not enough words in the bank: ROC_ID==" << info.roc
	      << ", BANK_ID==" << info.bankid <<" (0x" << std::hex << info.bankid << std::dec << ")"
	      << ", words_in_bank==" << info.length << std::endl;
    return kFALSE;
    
  } else {
    for (size_t i=0; i<info.length; i++){
      data = lbuff[i];
      if(decoder_index)             /* decoder type data word - set by decoder header word */
	{
	  type = 16;        /*  set decoder data words as type 16 */
	  new_type = 0;
	  if(decoder_index < num_decoder_words)
	    {
	      FillHDVariables(data, decoder_index - 1);
	      if(fDebug>1)
		printf("%8X - decoder data(%d) = %d\n", data, (decoder_index - 1),
		       data);
	      decoder_index++;
	    }
	  else                      /* last decoder word */

	    {
	      FillHDVariables(data, decoder_index - 1);
	      if(fDebug>1)
		printf("%8X - decoder data(%d) = %d\n", data, (decoder_index - 1),
		       data);
	      decoder_index = 0;
	      num_decoder_words = 1;
	    }
	}
      else                          /* normal typed word */
	{
	  if(data & 0x80000000)     /* data type defining word */
	    {
	      new_type = 1;
	      type = (data & 0x78000000) >> 27;
	    }
	  else                      /* data type continuation word */
	    {
	      new_type = 0;
	      type = type_last;
	    }
	  
	  switch (type)
	    {
	    case 0:         /* BLOCK HEADER */
	      slot_id_hd = (data & 0x7C00000) >> 22;
	      mod_id_hd = (data & 0x3C0000) >> 18;
	      n_evts = (data & 0x000FF);
	      blk_num = (data & 0x3FF00) >> 8;
	      if(fDebug>1)
		printf
		  ("%8X - BLOCK HEADER - slot = %d  id = %d  n_evts = %d  n_blk = %d\n",
		   data, slot_id_hd, mod_id_hd, n_evts,
		   blk_num);
	      break;
	      
	    case 1:         /* BLOCK TRAILER */
	      slot_id_tr = (data & 0x7C00000) >> 22;
	      n_words = (data & 0x3FFFFF);
	      if(fDebug>1)
		printf("%8X - BLOCK TRAILER - slot = %d   n_words = %d\n",
		       data, slot_id_tr, n_words);
	      break;
	      
	    case 2:         /* EVENT HEADER */
	      if(new_type)
		{
		  slot_id_ev_hd = (data & 0x07C00000) >> 22;
		  evt_num_1 = (data & 0x00000FFF);
		  trig_time = (data & 0x003FF000) >> 12;
		  if(fDebug>1)
		    printf
		      ("%8X - EVENT HEADER - slot = %d  evt_num = %d  trig_time = %d (%X)\n",
		       data, slot_id_ev_hd, evt_num_1, trig_time,
		       trig_time);
		}
	      break;
	      
	    case 3:         /* TRIGGER TIME */
	      
	      if(new_type)
		{
		  time_1 = (data & 0x7FFFFFF);
		  if(fDebug>1)
		    printf("%8X - TRIGGER TIME 1 - time = %X\n", data,
			   time_1);
		  time_now = 1;
		  time_last = 1;
		}
	      else
		{
		  if(time_last == 1)
		    {
		      time_2 = (data & 0xFFFFF);
		      if(fDebug>1)
			printf("%8X - TRIGGER TIME 2 - time = %X\n", data,
			       time_2);
		      time_now = 2;
		    }
		  else if(fDebug>1)
		    printf("%8X - TRIGGER TIME - (ERROR)\n", data);

		  time_last = time_now;
		}
          break;
	  
	    case 4:         /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;
	      
	    case 5:         /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;
	      
	    case 6:         /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;
	      
	    case 7:         /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;

	    case 8:         /* DECODER HEADER */
	      num_decoder_words = (data & 0x3F);    /* number of decoder words to follow */
	      num_words = num_decoder_words;
	      decoder_index = 1;    /* identify next word as a decoder data word */
	      if(fDebug>1)
		printf("%8X - DECODER HEADER = %d  (NUM DECODER WORDS = %d)\n",
		       data, type, num_words);
	      break;
	      
	    case 9:         /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;
	      
	    case 10:                /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;
	      
	    case 11:                /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;
	      
	    case 12:                /* UNDEFINED TYPE */
	      if(fDebug>1)
		printf("%8X - UNDEFINED TYPE = %d\n", data, type);
	      break;
	      
	    case 13:                /* END OF EVENT */
	      if(fDebug>1)
		printf("%8X - END OF EVENT = %d\n", data, type);
	      break;
	      
	    case 14:                /* DATA NOT VALID (no data available) */
	      slot_id_dnv = (data & 0x7C00000) >> 22;
	      if(fDebug>1)
		printf("%8X - DATA NOT VALID = %d  slot = %d\n", data,
		       type, slot_id_dnv);
	      break;
	      
	    case 15:                /* FILLER WORD */
	      slot_id_fill = (data & 0x7C00000) >> 22;
	      if(fDebug>1)
		printf("%8X - FILLER WORD = %d  slot = %d\n", data, type,
		       slot_id_fill);
	      break;
	    }
	  
	  type_last = type; /* save type of current data word */
	  
	}
    }
  }
  return kTRUE;
}

//____________________________________________________________________
void SBSHelicityDecoderBoard::FillHDVariables(uint32_t data,
					      uint32_t index){

  switch (index) {
  case 0:
    fSeed_Reported=data;
    break;
  case 1:
    fNum_TStable_Fall=data;
    break;
  case 2:
    fNum_TStable_Rise=data;
    break;
  case 3:
    fNum_Pattern_Sync=data;
    break;
  case 4:
    fNum_Pair_Sync=data;
    break;
  case 5:
    fTime_since_TStable=data;
    break;
  case 6:
    fTime_since_TSettle=data;
    break;
  case 7:
    fLast_Duration_TStable=data;
    break;
  case 8:
    fLast_Duration_TSettle=data;
    break;
  case 9:
    fPatternPhase       = (data>>8) & 0xff;
    fEventPolarity      = (data>>5) & 0x1;
    fReportedPatternHel = (data>>4) & 0x1;
    fBit_Helicity       = (data>>3) & 0x1;
    fBit_PairSync       = (data>>2) & 0x1;
    fBit_PatSync        = (data>>1) & 0x1;
    fBit_TStable        = data      & 0x1;
    break;
  case 10:
    fEvtHistory_PatSync=data;
    break;
  case 11:
    fEvtHistory_PairSync=data;
    break;
  case 12:
    fEvtHistory_ReportedHelicity=data;
    break;
  case 13:
    fPatHistory_ReportedHelicity=data;
    break;
  }
}


					      

//____________________________________________________________________
ClassImp(SBSHelicityDecoderBoard)

