///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSShowerBlock                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "SBSShowerBlock.h"

const double SBSShowerBlock::kBig = 1.e15;

SBSShowerBlock::SBSShowerBlock(Float_t x, Float_t y, 
                                   Float_t ped, Float_t gain, 
                                   Int_t row, Int_t col) 
{
    fX=x; fY=y; fPed=ped; fGain=gain; fRow=row; fCol=col;
}


SBSShowerBlock* SBSShowerBlock::operator=(SBSShowerBlock* rhs) {
    SetX(rhs->GetX());
    SetY(rhs->GetY());
    SetE(rhs->GetE());
    SetPed(rhs->GetPed());
    SetGain(rhs->GetGain());
    SetRow(rhs->GetRow());
    SetCol(rhs->GetCol());
    return this;
}

void SBSShowerBlock::ClearEvent() {
    fE = 0.; fStat = 0;
}



ClassImp(SBSShowerBlock)

