#ifndef ROOT_SBSScintBar
#define ROOT_SBSScintBar


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSScintBar                                                               //
//                                                                           //
// Class to represent a neutron bar                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class SBSScintBar;
class SBSScintPMT;

#include "TObject.h"
#include "SBSScintPMT.h"

class SBSScintBar : public TObject {

public: 
	SBSScintBar( Double_t x=0.,Double_t y=0.,Double_t z=0., 
		Double_t wx=1., Double_t wy=1., Double_t wz=1.,
		Double_t c=3.e8, Double_t att=0.,
		Double_t lgain=1., Int_t lped=0, Double_t lres=1.,
		Double_t loff=0., Double_t lwalk=0.,
		Double_t rgain=1., Int_t rped=0, Double_t rres=1.,
		Double_t roff=0., Double_t rwalk=0., Int_t barnum=0,
		Int_t llowlim=0, Int_t luplim=65535, Double_t lwrapa=0., 
		Int_t rlowlim=0, Int_t ruplim=65535, Double_t rwrapa=0. );

	virtual ~SBSScintBar();

	enum EBarType { kSpecial = -1, kUVA = 0, kGlasgow = 1, kCMU = 2 };
	void SetXPos( Double_t pos ) {fXPosPlane=pos;}
	void SetYPos( Double_t pos ) {fYPosPlane=pos;}
	void SetZPos( Double_t pos ) {fZPosPlane=pos;}
	void SetXWidth( Double_t w ) {fXWidth=w;}
	void SetYWidth( Double_t w ) {fYWidth=w;}
	void SetZWidth( Double_t w ) {fZWidth=w;}
	void SetC (Double_t c) {fc=c;}
	void SetAtt (Double_t a) {fatt=a;}
	void SetBarType( EBarType t) {fBarType=t;}
	void SetBarNum( Int_t barnum) {fBarNum=barnum;}
	void SetBarNum_nd( Int_t barnum_nd) {fBarNum_nd=barnum_nd;}

	Double_t GetXPos() const {return fXPosPlane;} 
	Double_t GetYPos() const {return fYPosPlane;} 
	Double_t GetZPos() const {return fZPosPlane;} 
	Double_t GetXWidth() const {return fXWidth;} 
	Double_t GetYWidth() const {return fYWidth;} 
	Double_t GetZWidth() const {return fZWidth;} 
	Double_t GetC() const {return fc;} 
	Double_t GetAtt() const {return fatt;}
	SBSScintPMT* GetLPMT() { return &fLPMT; }
	SBSScintPMT* GetRPMT() { return &fRPMT; }
	Double_t GetBarType() const {return fBarType;}
	Int_t GetBarNum() const {return fBarNum;}
	Int_t GetBarNum_nd() const {return fBarNum_nd;}

protected: 
	SBSScintPMT fLPMT;
	SBSScintPMT fRPMT;

	Double_t fXPosPlane;
	Double_t fYPosPlane;
	Double_t fZPosPlane;    // middle of the bar, relative to the center of the plane it belongs too

	Double_t fXWidth;  
	Double_t fYWidth;  
	Double_t fZWidth;       //full width (from end to end)

	Double_t fc;            //  effective speed of light  
	Double_t fatt;          //  attenuation
	EBarType fBarType;
	Int_t fBarNum;
	Int_t fBarNum_nd;

public:
	ClassDef(SBSScintBar,1) // Scintillator bar (PMT pointers/geometry) info.

};
////////////////////////////////////////////////////////////////////////////////

#endif
