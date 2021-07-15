#ifndef SBSSIMDATADECODER_H
#define SBSSIMDATADECODER_H

#include <vector>
#include <string>

#define SBS_MAX_ENCODER_WORDS 1024 // Max number of words the encoder can encode
#define SBS_NWORDS_MASK       0x3FFF
#define SBS_CHANNEL_MASK      0x3FFF // 0x1FF // 
#define SBS_TYPE_MASK         0xF // 0x1FF // 
#define SBS_CHANNEL_FIRST_BIT 14
#define SBS_TYPE_FIRST_BIT    28 // 23 // 
#define SBS_APV25_NCH 128 // Number of channels per APV25
#define SBS_MPD_NAPV25 15 // Number of AVP25's per MPD

//TODO: simplify this... I'm not convinced such a complexity is required anymore...

namespace SimEncoder {
  struct data {
    unsigned int channel;
  };

  struct adc_data : data {
    unsigned int integral;
  };
  
  //standard sample ADC structure
  struct sadc_data : adc_data {
    unsigned int integral;
    std::vector<unsigned int> samples;
  };

  /*
  struct fadc_data : adc_data {
    std::vector<unsigned int> samples;
  };
  */
  
  struct mpd_data : data {
    /*
    unsigned short nsamples; ///< Number of samples per channel
    unsigned short nstrips; ///< Number of strips in this block (128 for APV25)
    unsigned short mpd_id;
    unsigned short gem_id;
    unsigned short adc_id;
    unsigned short i2c;
    unsigned short pos;
    unsigned short invert;
    */
    std::vector<unsigned int> strips;
    std::vector<unsigned int> samples;
  };
  
  struct tdc_data : data {
    std::vector<unsigned int> time;
    unsigned int getTime(unsigned int t) { return time[t]&0x7FFFFFFF; }
    unsigned int getEdge(unsigned int t) { return (time[t]>>31)&0x1; }
    // Note: bit 31 will be the edge (1 for trail, 0 for lead)
  };
};

class SBSSimDataDecoder {
public:
  SBSSimDataDecoder(const char *enc_name, unsigned short enc_id);
  virtual ~SBSSimDataDecoder() {};
  
  /*
  // Encoders
  virtual bool EncodeADC(SimEncoder::adc_data data, unsigned int *enc_data,
      unsigned short &nwords) { return false; };
  virtual bool EncodeTDC(SimEncoder::tdc_data data, unsigned int *enc_data,
      unsigned short &nwords) { return false; };
  virtual bool EncodeFADC(SimEncoder::fadc_data data, unsigned int *enc_data,
      unsigned short &nwords) { return false; };
  virtual bool EncodeMPD(SimEncoder::mpd_data data, unsigned int *enc_data,
      unsigned short &nwords) { return false; };
  */
  // Decoders
  virtual bool DecodeADC(SimEncoder::adc_data &data,
			 const unsigned int *enc_data,
			 unsigned short nwords) { return false; }
  virtual bool DecodeTDC(SimEncoder::tdc_data &data,
			 const unsigned int *enc_data,
			 unsigned short nwords) { return false; };
  virtual bool DecodeSADC(SimEncoder::sadc_data &data,
			  const unsigned int *enc_data,
			  unsigned short nwords) { return false; }
  virtual bool DecodeMPD(SimEncoder::mpd_data &data,
			 const unsigned int *enc_data,
			 unsigned short nwords) { return false; }
  /*
  virtual bool DecodeFADC(SimEncoder::sadc_data &data,
      const unsigned int *enc_data,unsigned short nwords) { return false; }
  virtual bool DecodeMPD(SimEncoder::sadc_data &data,
      const unsigned int *enc_data,unsigned short nwords) { return false; }
  virtual bool DecodeFADC(SimEncoder::fadc_data &data,
      const unsigned int *enc_data,unsigned short nwords) { return false; }
  virtual bool DecodeMPD(SimEncoder::mpd_data &data,
      const unsigned int *enc_data,unsigned short nwords) { return false; }
  */
  // Capabilities
  virtual bool IsADC() { return false; }
  virtual bool IsTDC() { return false; }
  virtual bool IsSADC() { return false; }
  //virtual bool IsFADC() { return false; }
  virtual bool IsMPD() { return false; }

  unsigned short GetId() { return fEncId; }
  std::string GetName() { return fName; }

  static SBSSimDataDecoder* GetEncoderByName(const char *enc_name);
  static SBSSimDataDecoder* GetEncoder(unsigned short id);
  static unsigned int MakeBitMask(unsigned short bits);

  static unsigned int EncodeHeader(unsigned short type, unsigned short mult,
      unsigned int nwords);
  static void DecodeHeader(unsigned int hdr, unsigned short &type, unsigned short &ch,
      unsigned int &nwords);
  static unsigned short DecodeChannel(unsigned int hdr);
  static unsigned short DecodeType(unsigned int hdr);
  static unsigned short DecodeNwords(unsigned int hdr);

protected:
  std::string fName;
  unsigned short fEncId;

private:
  static std::vector<SBSSimDataDecoder*> fEncoders;
};

// Generic TDC encoder
class SBSSimTDCEncoder : public SBSSimDataDecoder {
public:
  SBSSimTDCEncoder(const char *enc_name, unsigned short enc_id,
      unsigned short bits, unsigned short edge_bit);
  virtual ~SBSSimTDCEncoder() {};

  // Overloaded functions
  //virtual bool EncodeTDC(SimEncoder::tdc_data data, unsigned int *enc_data,
  //  unsigned short &nwords);
  virtual bool DecodeTDC(SimEncoder::tdc_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsTDC() { return true; }

protected:
  unsigned short fBits;
  unsigned short fEdgeBit;
  unsigned short fBitMask;
};

// Generic ADC encoder
class SBSSimADCEncoder : public SBSSimDataDecoder {
public:
  SBSSimADCEncoder(const char *enc_name, unsigned short enc_id,
      unsigned short bits);
  virtual ~SBSSimADCEncoder() {};

  //virtual bool EncodeADC(SimEncoder::adc_data data, unsigned int *enc_data,
  //  unsigned short &nwords);
  virtual bool DecodeADC(SimEncoder::adc_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsADC() { return true; }

protected:
  unsigned short fBits;
  unsigned short fBitMask;
};

//generic sample ADC encoder, meant to replace FADC250 and MPD
class SBSSimSADCEncoder : public SBSSimADCEncoder {
public:
  SBSSimSADCEncoder(const char *enc_name, unsigned short enc_id);
  virtual ~SBSSimSADCEncoder() {};

  virtual bool DecodeSADC(SimEncoder::sadc_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsSADC() { return true; }
};

/*
// JLab FADC 250 in multi sample ADC mode
class SBSSimFADC250Encoder : public SBSSimADCEncoder {
public:
  SBSSimFADC250Encoder(const char *enc_name, unsigned short enc_id);
  virtual ~SBSSimFADC250Encoder() {};

  //virtual bool EncodeFADC(SimEncoder::fadc_data data, unsigned int *enc_data,
  //  unsigned short &nwords);
  //virtual bool DecodeFADC(SimEncoder::fadc_data &data,
  virtual bool DecodeFADC(SimEncoder::sadc_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsFADC() { return true; }

  //private:
  //unsigned int EncodeSingleSample(unsigned int dat);
  //void UnpackSamples(unsigned int enc_data,unsigned int *buff,
  //  bool *overflow, bool *valid);
};
*/

// MPD
class SBSSimMPDEncoder : public SBSSimDataDecoder {
public:
  SBSSimMPDEncoder(const char *enc_name, unsigned short enc_id);
  virtual ~SBSSimMPDEncoder() {};

  //virtual bool EncodeMPD(SimEncoder::mpd_data data, unsigned int *enc_data,
  //  unsigned short &nwords);
  virtual bool DecodeMPD(SimEncoder::mpd_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsMPD() { return true; }

  //void EncodeMPDHeader(SimEncoder::mpd_data data, unsigned int *enc_data,
  //  unsigned short &nwords);
  //void DecodeMPDHeader(const unsigned int *hdr, SimEncoder::mpd_data &data);
  /*
protected:
  unsigned int fChannelBitMask;
  unsigned int fDataBitMask;
  unsigned int fOverflowBitMask;
  unsigned int fSampleBitMask;
  unsigned int fValidBitMask;
  */
  //protected:
  //unsigned int EncodeSingleSample(unsigned int dat);
  //void UnpackSamples(unsigned int enc_data,unsigned int *buff,
  //  bool *overflow, bool *valid);

};
/**/

#endif // SBSSIMDATAENCODER_H
