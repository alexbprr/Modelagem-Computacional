#ifndef H_MULTIARRAY
#define H_MULTIARRAY
#include "Array.h"

/*
   The MultiArray class is an abstraction for a multi dimensional array that is implemented as a 
   contiguos one dimensional array. This class inherits the Array class and defines an attribute 
   to store the size of each spatial dimension.
   An one dimensional array can be implemented by defining the size of the second and third dimensions to be zero.
   A two dimensional array can be implemented by defining the size of the third dimension to be zero.
*/
class MultiArray : public virtual Array {
    int dl[3]; // length of each dimension
public:
    MultiArray(int l[3]): Array(l[0]*l[1]*l[2]){
        dl[0] = l[0], dl[1] = l[1], dl[2] = l[2];        
    }
    ~MultiArray(){}
    int index (int i, int j, int k) const { 
        return k * dl[0] * dl[1] + j * dl[0] + i;
    }    
    double operator () (int i, int j, int k) const { return v[index(i,j,k)]; }
    double &operator () (int i, int j, int k){ return v[index(i,j,k)]; }
};

#endif 