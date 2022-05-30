#ifndef H_MULTIVALUEARRAY
#define H_MULTIVALUEARRAY
#include <stdlib.h>

class MultiValueArray {
protected:
    double** v;
    int ndim;
    int length;

public:
    MultiValueArray(int n, int l){ 
        ndim = n;
        length = l;         
        v = (double**) calloc(n, sizeof(double*));
        for (int i = 0; i < n; i++)
            v[i] = (double*) calloc(l, sizeof(double));
    }
    virtual ~MultiValueArray(){
        for (int i = 0; i < ndim; i++)
            delete v[i];
        delete v; 
    }
    double operator () (int i, int j) const{ return v[i][j]; }
    double &operator () (int i, int j) { return v[i][j]; }
    int getLength(){ return length; }
    int getNumberOfDims(){ return ndim; }
    void setValue(int i, int j, double va){ v[i][j] = va; }
};

#endif 