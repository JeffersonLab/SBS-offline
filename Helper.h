#ifndef SBS_Offline_Helper_h
#define SBS_Offline_Helper_h

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Helper                                                                    //
//                                                                           //
// Helper classes and functions                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>
#include <algorithm>

//___________________________________________________________________________
struct DeleteObject {
  template<typename T>
  void operator()( const T* ptr ) const { delete ptr; }
};

//___________________________________________________________________________
template<typename Container>
inline void DeleteContainer( Container& c ) {
  // Delete all elements of given container of pointers
  for_each(c.begin(), c.end(), DeleteObject());
  c.clear();
}

//___________________________________________________________________________
template<typename ContainerOfContainers>
inline void DeleteContainerOfContainers( ContainerOfContainers& cc ) {
  // Delete all elements of given container of containers of pointers
  for_each(cc.begin(), cc.end(),
           DeleteContainer<typename ContainerOfContainers::value_type>);
  cc.clear();
}

//___________________________________________________________________________
//inline Int_t NumberOfSetBits( UInt_t v ) {
//  // Count number of bits set in 32-bit integer. From
//  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
//
//  v = v - ((v >> 1) & 0x55555555);
//  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
//  v = (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
//  return static_cast<Int_t>(v);
//}

//___________________________________________________________________________
// Simple and quick routine to init and clear most vectors
// (of integers, doubles, doubles, etc...)
// Reset/Init 1D vector
template<class T>
void InitVector(std::vector<T> &vec, T val = 0, size_t n = 0) {
  vec.resize(n);
  ResetVector(vec,val);
}

//___________________________________________________________________________
template<class T>
void ResetVector(std::vector<T> &vec, T val = 0, size_t n = 0) {
  if(n > 0) {
    vec.clear();
    vec.resize(n);
  }
  vec.assign(vec.size(),val);
}

//___________________________________________________________________________
// Reset 2D vector
template<class T>
void ResetVector(std::vector<std::vector<T> > &vec, T val = 0,
                 size_t nr = 0, size_t nc = 0) {
  if(nr > 0) {
    vec.clear();
    vec.resize(nr);
  }
  for(size_t i = 0; i < vec.size(); i++) {
    ResetVector(vec[i],val,nc);
  }
}

///////////////////////////////////////////////////////////////////////////////

#endif
