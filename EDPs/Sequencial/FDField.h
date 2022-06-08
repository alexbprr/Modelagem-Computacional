#ifndef H_FDFIELD
#define H_FDFIELD
#include <string> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <unordered_map>
#include <vector> 
#include <math.h>
using namespace std; 

#include "MultiArray.h"
#include "IntMultiArray.h"

enum Dim { d1D, d2D, d3D};

class FDField {
protected:
    std::string name; 
    std::string filename; 
    int dim[3]; //Size of x, y and z dimensions
    enum Dim Ndim; //Number of dimensions (1D, 2D or 3D)
    double dx, dy, dz; //Spatial discretization for each dimension 
    int Nx,Ny,Nz; //Number of points in x, y and z dimensions 
    Array *xarray; //Array that stores the values in the x dimension 
    Array *yarray; //Array that stores the values in the y dimension 
    Array *zarray; //Array that stores the values in the z dimension 
    MultiArray *values; //Pointer to array that stores the field values 
    MultiArray *nextvalues;        
    IntMultiArray *sources; //sources(i,j,k) == 1: there is a source in index i,j,k
    IntMultiArray *sinks; //sink(i,j,k) == 1: there is a sink in index i,j,k

    double nonlinearDiffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, double D, string W, int x, int y, int z);
    double classicDiffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, double D, int x, int y, int z);
    virtual double diffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z); 
    double classicChemotaxis(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, std::string u, std::string c, double X, int x, int y, int z);
    virtual double chemotaxis(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z);
    virtual double reactions(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z);        
    virtual bool update(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z, double dt);     

public:    

    FDField(string n, int d[3], double deltas[3]); 
    virtual ~FDField();     
    double operator () (int i, int j, int k) { return (*values)(i,j,k); }    
    double getValue(double x, double y, double z);
    void setValue (int i, int j, int k, double va) { (*nextvalues)(i,j,k) = va; }
    MultiArray* getValues(){ return values;}  
    std::string getFileName(){ return filename; }     
    virtual bool checkErrors(double v);
    bool convertPosToIndex(double x, double y, double z, int *index);
    void setInitialValue(double x, double y, double z, double v);
    void addSource(double x, double y, double z);
    int source(double x, double y, double z); 
    int hasSource(int x, int y, int z);  
    void addSink(double x, double y, double z);    
    int sink(double x, double y, double z);    
    void printSources();
    void printSinks();
    virtual void update(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, double dt); 
    void saveX(){
        std::ofstream fp;
        std::stringstream ss;
        fp.open("x.csv");
        ss << "x" << endl;
        for (int i = 0; i < Nx; i++)
            ss << (*xarray)(i) << endl;
        fp << ss.str();
        fp.close();
    }
    void saveY(){
        std::ofstream fp;
        std::stringstream ss;
        fp.open("y.csv");
        ss << "y" << endl;
        for (int i = 0; i < Ny; i++)
            ss << (*yarray)(i) << endl;        
        fp << ss.str();
        fp.close();
    }
    void saveZ(){
        std::ofstream fp;
        std::stringstream ss;
        fp.open("z.csv");
        ss << "z" << endl;
        for (int i = 0; i < Nz; i++)
            ss << (*zarray)(i) << endl;
        fp << ss.str();
        fp.close();
    }
    void saveXYZ(){
        saveX();
        saveY();
        saveZ();
    }
    void save(double t, std::string fpath);     
};

#endif 