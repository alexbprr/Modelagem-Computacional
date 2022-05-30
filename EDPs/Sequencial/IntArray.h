#ifndef H_INTARRAY
#define H_INTARRAY
#include <stdlib.h>

class IntArray {
protected:
    int* v;
    int length;
public:
    IntArray(int l){ 
        length = l; 
        v = (int*) calloc(l, sizeof(int));         
    }
    virtual ~IntArray(){
        delete v; 
    }
    int operator () (int i) const{ return v[i]; }
    int &operator () (int i) { return v[i]; }
    int getLength(){ return length; }
    void setValue(int i, int va){ v[i] = va; }
};

#endif 