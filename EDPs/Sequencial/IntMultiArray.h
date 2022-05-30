#ifndef H_INTMULTIARRAY
#define H_INTMULTIARRAY
#include "IntArray.h"

class IntMultiArray : public virtual IntArray {
    int dl[3]; // length of each dimension
public:
    IntMultiArray(int l[3]): IntArray(l[0]*l[1]*l[2]){
        dl[0] = l[0], dl[1] = l[1], dl[2] = l[2];        
    }
    ~IntMultiArray(){}
    int index (int i, int j, int k) const { 
        return k * dl[0] * dl[1] + j * dl[0] + i;
    }    
    int operator () (int i, int j, int k) const { return v[index(i,j,k)]; }
    int &operator () (int i, int j, int k){ return v[index(i,j,k)]; }
};

#endif 