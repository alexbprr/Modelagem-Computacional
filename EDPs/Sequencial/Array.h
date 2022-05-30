#ifndef H_ARRAY
#define H_ARRAY
#include <stdlib.h>

class Array {
protected:
    double* v;
    int length;
public:
    Array(int l){ 
        length = l; 
        v = (double*) calloc(l, sizeof(double));         
    }
    virtual ~Array(){
        delete v; 
    }
    double operator () (int i) const{ return v[i]; }
    double &operator () (int i) { return v[i]; }
    int getLength(){ return length; }
    void setValue(int i, double va){ v[i] = va; }
};

#endif 