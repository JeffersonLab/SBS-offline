#include "SBSSimDataDecoder.h"
#include <TString.h>
#include <iostream>

#define SBS_MPD_NSAMPLES_BIT  27
#define SBS_MPD_NSTRIPS_BIT   15
#define SBS_MPD_POS_BIT       7
#define SBS_MPD_INV_BIT       6
#define SBS_MPD_I2C_BIT       0

#define SBS_MPD_ADC_ID_BIT    20
#define SBS_MPD_MPD_ID_BIT    10
#define SBS_MPD_GEM_ID_BIT    0

// Word 1
#define SBS_MPD_NSAMPLES_MASK 0xF8000000
#define SBS_MPD_NSTRIPS_MASK  0x07FF8000
#define SBS_MPD_POS_MASK      0x00007F80
#define SBS_MPD_INV_MASK      0x10000040
#define SBS_MPD_I2C_MASK      0x0000003F
// Word 2
#define SBS_MPD_ADC_ID_MASK   0x3FF00000
#define SBS_MPD_MPD_ID_MASK   0x000FFC00
#define SBS_MPD_GEM_ID_MASK   0x000003FF

// This one is static, so define it again here
std::vector<SBSSimDataDecoder*> SBSSimDataDecoder::fEncoders;

SBSSimDataDecoder* SBSSimDataDecoder::GetEncoderByName(
    const char *enc_name)
{
  if(fEncoders.empty()) { // First generate the list of known encoders!!
    unsigned short ids = 1;
    // TDCs
    fEncoders.push_back(new SBSSimTDCEncoder("caen1190",ids++,19,26));
    fEncoders.push_back(new SBSSimTDCEncoder("lecroy1877",ids++,16,16));
    fEncoders.push_back(new SBSSimTDCEncoder("vetroc",ids++,16,26));
    fEncoders.push_back(new SBSSimTDCEncoder("f1tdc",ids++,16,31));
    // ADCs
    //fEncoders.push_back(new SBSSimFADC250Encoder("fadc250",ids++));
    fEncoders.push_back(new SBSSimSADCEncoder("fadc250",ids++));
    fEncoders.push_back(new SBSSimADCEncoder("adc",ids++,12));
    fEncoders.push_back(new SBSSimADCEncoder("lecroy1881",ids++,14));
    fEncoders.push_back(new SBSSimADCEncoder("caen792",ids++,12));
    //fEncoders.push_back(new SBSSimSADCEncoder("mpd",ids++));
    fEncoders.push_back(new SBSSimMPDEncoder("mpd",ids++));
  }

  TString name(enc_name);
  for(std::vector<SBSSimDataDecoder*>::iterator it = fEncoders.begin();
      it != fEncoders.end(); it++) {
    if(name.CompareTo((*it)->GetName(),TString::kIgnoreCase)==0)
      return *it;
  }

  return 0;
}

SBSSimDataDecoder* SBSSimDataDecoder::GetEncoder(unsigned short id)
{
  for(std::vector<SBSSimDataDecoder*>::iterator it = fEncoders.begin();
      it != fEncoders.end(); it++) {
    if((*it)->GetId() == id)
      return *it;
  }
  return 0;
}

SBSSimDataDecoder::SBSSimDataDecoder(const char *enc_name,
    unsigned short enc_id) : fName(enc_name), fEncId(enc_id)
{
}

SBSSimTDCEncoder::SBSSimTDCEncoder(const char *enc_name,
    unsigned short enc_id, unsigned short bits, unsigned short edge_bit)
  : SBSSimDataDecoder(enc_name,enc_id), fBits(bits), fEdgeBit(edge_bit)
{
  fBitMask = MakeBitMask(fBits);
}

SBSSimADCEncoder::SBSSimADCEncoder(const char *enc_name,
    unsigned short enc_id, unsigned short bits)
  : SBSSimDataDecoder(enc_name,enc_id), fBits(bits)
{
  fBitMask = MakeBitMask(fBits);
}

SBSSimSADCEncoder::SBSSimSADCEncoder(const char *enc_name,
    unsigned short enc_id) : SBSSimADCEncoder(enc_name,enc_id,12)
{
}

/*
SBSSimFADC250Encoder::SBSSimFADC250Encoder(const char *enc_name,
    unsigned short enc_id) : SBSSimADCEncoder(enc_name,enc_id,12)
{
}
*/

unsigned int SBSSimDataDecoder::MakeBitMask(unsigned short bits)
{
  unsigned int mask = 0;
  for(unsigned short b = 0; b < bits; b++) {
    mask |= 1<<b;
  }
  return mask;
}

/*
bool SBSSimADCEncoder::EncodeADC(SimEncoder::adc_data data,
    unsigned int *enc_data,unsigned short &nwords)
{
  nwords = 0;
  enc_data[nwords++] = data.integral&fBitMask;
  return (nwords>0);
}
*/
/*
bool SBSSimTDCEncoder::EncodeTDC(SimEncoder::tdc_data data,
    unsigned int *enc_data,unsigned short &nwords)
{
  // Generic TDC encoder where the lowest n-bits are the time
  // measurement, and a single edge bit specifies either 0: lead, 1: trail
  // (We ignore the channel bit because it goes unused in this simulation)
  nwords = 0;
  for(std::vector<unsigned int>::iterator it = data.time.begin();
      it != data.time.end(); it++) {
    enc_data[nwords++] = ((*it)&fBitMask) | ((((*it)>>31)&0x1)<<fEdgeBit);
  }
  return (nwords>0);
}
*/
/*
bool SBSSimFADC250Encoder::EncodeFADC(SimEncoder::fadc_data data,
    unsigned int *enc_data, unsigned short &nwords)
{
  // Word 1:
  //   (31) = 1 (<-- Skipping, not useful for Sim data)
  //   (30 – 27) = 4 (<-- I'm going to skip this, not sure it is useful)
  //   (26 – 23) = channel number (0 – 15) (<-- Also skipping)
  //   (22 – 12) = reserved (read as 0)
  //   (11 – 0) = window width (in number of samples)
  // Words 2 - N:
  //   (31) = 0
  //   (30) = reserved (read as 0)
  //   (29) = sample x not valid
  //   (28 – 16) = ADC sample x (includes overflow bit)
  //   (15 – 14) = reserved (read as 0)
  //   (13) = sample x + 1 not valid
  //   (12 – 0) = ADC sample x + 1 (includes overflow bit)
  nwords = 0;
  unsigned int nsamps = data.samples.size() >= 0xFFF ? 0xFFF :
    data.samples.size() ;
  enc_data[nwords++] = nsamps;
  unsigned int s = 0;
  unsigned int buff[2] = {0,0};
  for(s = 0; s < nsamps-1; s+=2) {
    buff[0] = EncodeSingleSample(data.samples[s]);
    buff[1] = EncodeSingleSample(data.samples[s+1]);
    enc_data[nwords++] = (buff[0]<<16) | buff[1];
  }
  if( s < nsamps ) { // Still have one more sample to process
    buff[0] = EncodeSingleSample(data.samples[s]);
    buff[1] = 0x2000; // Mark last sample in this two-sample word as not valid
    enc_data[nwords++] = (buff[0]<<16) | buff[1];
  }
  return (nwords>1);
}
*/

bool SBSSimADCEncoder::DecodeADC(SimEncoder::adc_data &data,
      const unsigned int *enc_data,unsigned short nwords)
{
  if(nwords>1)
    return false;
  unsigned short nread = 0;
  data.integral = enc_data[nread++]&fBitMask;
  //std::cout << enc_data[0] << " " << data.integral << std::endl; 
  return nread==nwords;
}

bool SBSSimTDCEncoder::DecodeTDC(SimEncoder::tdc_data &data,
				 const unsigned int *enc_data,
				 unsigned short nwords)
{
  for(unsigned short n = 0; n < nwords; n++) {
    //std::cout << "n = " << n << ": encoded data " << enc_data[n] << " edge bit " << fEdgeBit << std::endl;
    //data.time.push_back(((enc_data[n]>>fEdgeBit)<<31) |
    //  (enc_data[n]&fBitMask));
    data.time.push_back(enc_data[n]);
    //std::cout << " decoded data " << (((enc_data[n]>>fEdgeBit)<<31) | (enc_data[n]&fBitMask)) << " edge bit " << fEdgeBit << std::endl;
    //std::cout << data.time[n] << " " << data.getTime(n) << " " << data.getEdge(n) << std::endl;
  }
  return !data.time.empty();
}

bool SBSSimSADCEncoder::DecodeSADC(SimEncoder::sadc_data &data,
    const unsigned int *enc_data,unsigned short nwords)
{
  data.integral = 0;
  data.samples.clear(); // Clear out any samples already in the data
  //std::cout << " nwords: " << nwords << " => ";
  for(unsigned short i = 0; i<nwords; i++){
    data.integral+=enc_data[i];
    data.samples.push_back(enc_data[i]);
    //std::cout << enc_data[i] << " ";
  }
  //std::cout << std::endl;
  return true;
}

/*
bool SBSSimFADC250Encoder::DecodeFADC(//SimEncoder::fadc_data &data,
				      SimEncoder::sadc_data &data,
    const unsigned int *enc_data,unsigned short nwords)
{
  for(unsigned short i = 0; i<nwords; i++){
    data.samples.push_back(enc_data[i]);
  }
  
  //OK, so the stuff below is flat-out out of date with the new digitization paradigm

  int nsamples = enc_data[0]&0xFFF;
  int nsamples_read = 0;

  std::cout << enc_data[0] << " " << nsamples << std::endl;

  unsigned int buff[2] = {0,0};
  bool overflow[2] = { false, false};
  bool valid[2] = {false, false};
  for(unsigned short n = 1; n < nwords; n++) {
    UnpackSamples(enc_data[n],buff,overflow,valid);
    for(short k = 0; k < 2; k++) {
      if(valid[k]) {
        data.samples.push_back(buff[k]);
        nsamples_read++;
      }
    }
  }

  if(nsamples_read != nsamples) {
    std::cerr << "DecodeFADC: Error, number of samples read (" << nsamples_read
      << "), does not match number of expected samples (" << nsamples
      << ")." << std::endl;
    return false;
  }

  return true;

}
*/
/*
unsigned int SBSSimFADC250Encoder::EncodeSingleSample(unsigned int dat)
{
  if(dat&0xFFFFF000) { // Data too large, turn on overflow
    dat = 0x1FFF;
  }
  return dat&0x1FFF;
}
*/
/*
void SBSSimFADC250Encoder::UnpackSamples(unsigned int enc_data,
    unsigned int *buff, bool *overflow, bool *valid)
{
  unsigned int tmp;
  for(int k = 0; k < 2; k++) {
    tmp = (k==0 ? (enc_data>>16) : enc_data)&0x3FFF;
    buff[k] = tmp&0xFFF;
    overflow[k] = tmp&0x1000;
    valid[k] = !(tmp&0x2000);
  }
}
*/

unsigned int SBSSimDataDecoder::EncodeHeader(unsigned short type,
    unsigned short mult, unsigned int nwords)
{
  // First word bits
  // 31-26: encoder type
  // 25-14: channel multiplier (to be converted to local channel by SimDecoder)
  // 13-0 : number of words that follow
  return ((type&SBS_TYPE_MASK)<<SBS_TYPE_FIRST_BIT) |
    ((mult&SBS_CHANNEL_MASK)<<SBS_CHANNEL_FIRST_BIT) |
    (nwords&SBS_NWORDS_MASK);
}

void SBSSimDataDecoder::DecodeHeader(unsigned int hdr, unsigned short &type, unsigned short &ch,
    unsigned int &nwords)
{
  type = DecodeType(hdr);
  ch = DecodeChannel(hdr);
  nwords = DecodeNwords(hdr);
}

unsigned short SBSSimDataDecoder::DecodeChannel(unsigned int hdr) {
  return (hdr>>SBS_CHANNEL_FIRST_BIT)&SBS_CHANNEL_MASK;
}

unsigned short SBSSimDataDecoder::DecodeType(unsigned int hdr) {
  return (hdr>>SBS_TYPE_FIRST_BIT)&SBS_TYPE_MASK;
}

unsigned short SBSSimDataDecoder::DecodeNwords(unsigned int hdr) {
  return hdr&SBS_NWORDS_MASK;
}

SBSSimMPDEncoder::SBSSimMPDEncoder(const char *enc_name,
    unsigned short enc_id) : SBSSimDataDecoder(enc_name,enc_id)
{
/*
  fChannelBitMask  = MakeBitMask(8);
  fDataBitMask     = MakeBitMask(12);
  fOverflowBitMask = (1<<12);
  fSampleBitMask   = fDataBitMask|fOverflowBitMask;
  fValidBitMask    = (1<<13);
*/
}
/*
bool SBSSimMPDEncoder::EncodeMPD(SimEncoder::mpd_data data,
    unsigned int *enc_data, unsigned short &nwords)
{
  // Word 1:
  //   31-27 = nsamples per strip (i.e. 3 or 6)
  //   26-15 = nstrips
  //   14-7  = pos
  //    6    = inv
  //    5-0  = i2c
  // Word 2:
  //   31-30 = not used
  //   29-20 = adc_id
  //   19-10 = mpd_id
  //    9-0  = gem_id
  // Word 3 - N: (Same logic as FADC250)
  //   31-30 = Not used for simulation
  //   29    = sample x not valid
  //   28–16 = ADC sample x (includes overflow bit)
  //   15–14 = Not used for simulation
  //   13 = sample x + 1 not valid
  //   12–0 = ADC sample x + 1 (includes overflow bit)
  nwords = 0;
  unsigned int nsamps = data.nsamples*data.nstrips;
  if(nsamps != data.samples.size())
    return false;
  EncodeMPDHeader(data,enc_data,nwords);
  unsigned int s = 0;
  unsigned int buff[2] = {0,0};
  for(s = 0; s < nsamps-1; s+=2) {
    buff[0] = EncodeSingleSample(data.samples[s]);
    buff[1] = EncodeSingleSample(data.samples[s+1]);
    enc_data[nwords++] = (buff[0]<<16) | buff[1];
  }
  if( s < nsamps ) { // Still have one more sample to process
    buff[0] = EncodeSingleSample(data.samples[s]);
    buff[1] = fValidBitMask; // Mark last sample in this two-sample word as not valid
    enc_data[nwords++] = (buff[0]<<16) | buff[1];
  }
  return (nwords>1);
}
*/
/*
unsigned int SBSSimMPDEncoder::EncodeSingleSample(unsigned int dat)
{
  if(dat>fDataBitMask) { // Data too large, turn on overflow
    dat |= fOverflowBitMask;
  }
  return dat&fSampleBitMask;
}
*/
/*
void SBSSimMPDEncoder::UnpackSamples(unsigned int enc_data,
    unsigned int *buff, bool *overflow, bool *valid)
{
  unsigned int tmp;
  for(int k = 0; k < 2; k++) {
    tmp = k==0 ? (enc_data>>16) : enc_data;
    buff[k] = tmp&fDataBitMask;
    overflow[k] = tmp&fOverflowBitMask;
    valid[k] = !(tmp&fValidBitMask);
  }
}
*/
/*
void SBSSimMPDEncoder::EncodeMPDHeader(SimEncoder::mpd_data data,
    unsigned int *enc_data, unsigned short &nwords)
{
  enc_data[nwords++] =
    ((data.nsamples<<SBS_MPD_NSAMPLES_BIT)&SBS_MPD_NSAMPLES_MASK) |
    ((data.nstrips<<SBS_MPD_NSTRIPS_BIT)&SBS_MPD_NSTRIPS_MASK) |
    ((data.pos<<SBS_MPD_POS_BIT)&SBS_MPD_POS_MASK) |
    ((data.invert<<SBS_MPD_INV_BIT)&SBS_MPD_INV_MASK) |
    (data.i2c&SBS_MPD_I2C_MASK);
  enc_data[nwords++] =
    ((data.adc_id<<SBS_MPD_ADC_ID_BIT)&SBS_MPD_ADC_ID_MASK) |
    ((data.mpd_id<<SBS_MPD_MPD_ID_BIT)&SBS_MPD_MPD_ID_MASK) |
    (data.gem_id&SBS_MPD_GEM_ID_BIT);
}
*/
/*
void SBSSimMPDEncoder::DecodeMPDHeader(const unsigned int *hdr,
    SimEncoder::mpd_data &data)
{
  // Word 1
  data.nsamples = (*hdr&SBS_MPD_NSAMPLES_MASK)>>SBS_MPD_NSAMPLES_BIT;
  data.nstrips  = (*hdr&SBS_MPD_NSTRIPS_MASK)>>SBS_MPD_NSTRIPS_BIT;
  data.pos      = (*hdr&SBS_MPD_POS_MASK)>>SBS_MPD_POS_BIT;
  data.invert   = (*hdr&SBS_MPD_INV_MASK)>>SBS_MPD_INV_BIT;
  data.i2c      = (*hdr++&SBS_MPD_I2C_MASK);
  // Word 2
  data.adc_id   = (*hdr&SBS_MPD_ADC_ID_MASK)>>SBS_MPD_ADC_ID_BIT;
  data.mpd_id   = (*hdr&SBS_MPD_MPD_ID_MASK)>>SBS_MPD_MPD_ID_BIT;
  data.gem_id   = (*hdr&SBS_MPD_GEM_ID_MASK)>>SBS_MPD_GEM_ID_BIT;
}
*/
bool SBSSimMPDEncoder::DecodeMPD(SimEncoder::mpd_data &data,
				 const unsigned int *enc_data,unsigned short nwords)
{
  data.samples.clear(); // Clear out any samples already in the data
  //std::cout << "nwords " << nwords << std::endl;
  //data.nstrips = nwords/data.nsamples;
  for(unsigned short i = 0; i<nwords; i++){
    //std::cout << enc_data[i] << " ";
    
    data.strips.push_back(enc_data[i]/8192);
    data.samples.push_back(enc_data[i]%8192);
  }
  //std::cout << std::endl;
  /*
  if(nwords<=1) {
    std::cerr << "Error, not enough words to read. Expected more than one,"
      << " got only: " << nwords << std::endl;
    return false;
  }

  // First, decode the header
  DecodeMPDHeader(enc_data,data);
  int nsamples_read = 0;
  data.samples.clear(); // Clear out any samples already in the data

  unsigned int buff[2] = {0,0};
  bool overflow[2] = { false, false};
  bool valid[2] = {false, false};
  for(unsigned short n = 2; n < nwords; n++) {
    UnpackSamples(enc_data[n],buff,overflow,valid);
    for(short k = 0; k < 2; k++) {
      if(valid[k]) {
        data.samples.push_back(buff[k]);
        nsamples_read++;
      }
    }
  }

  if(nsamples_read != data.nstrips*data.nsamples) {
    std::cerr << "Error, number of samples read (" << std::dec << nsamples_read
      << "), does not match number of expected samples ("
      << data.nstrips*data.nsamples << ")." << std::endl;
    return false;
  }
  */
  return true;

}
/*
*/
