#ifndef ROOT_SBSShowerBlock
#define ROOT_SBSShowerBlock

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// SBSShowerBlock                                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class SBSShowerBlock : public TObject {

public:

    SBSShowerBlock() {}
    SBSShowerBlock(Float_t x, Float_t y, Float_t ped,
        Float_t gain, Int_t row, Int_t col);

    Float_t GetX()    const {return fX;   }
    Float_t GetY()    const {return fY;   }
    Float_t GetE()    const {return fE;   }
    Float_t GetPed()  const {return fPed; }
    Float_t GetGain() const {return fGain;}
    Int_t   GetRow()  const {return fRow; }
    Int_t   GetCol()  const {return fCol; }
    Int_t   GetStat() const {return fStat;}

    void SetX(Float_t var)    { fX = var;    }
    void SetY(Float_t var)    { fY = var;    }
    void SetE(Float_t var)    { fE = var;    }
    void SetPed(Float_t var)  { fPed = var;  }
    void SetGain(Float_t var) { fGain = var; }
    void SetRow(Int_t var)    { fRow = var;  }
    void SetCol(Int_t var)    { fCol = var;  }
    void SetStat(Int_t var)   { fStat = var; }

    SBSShowerBlock* operator=(SBSShowerBlock* rhs);

    void ClearEvent();

private:

    static const double kBig; //  = 1e15;

    Float_t fX;       // x position of the center
    Float_t fY;       // y position of the center
    Float_t fPed;     // Pedestal
    Float_t fGain;    // Gain

    Float_t fE;       // Energy deposit in block

    Int_t   fRow;     // Row of the block
    Int_t   fCol;     // Column of the block

    Int_t   fStat;    // Status: 0: not seen, 1: seen, 2: local max

    ClassDef(SBSShowerBlock,0)   // Generic shower block class
};

#endif
