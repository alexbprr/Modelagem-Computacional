#include "PdeSolver.h"
#include <random>

/*Use (*xarray)(x),(*yarray)(y) and (*xarray)(z) to access the values of the coordinates of a point
*/ 
//Create a field for each population of the system
double betap = 1,
     mp = 0.2, 
     Dp = 0.1, 
     Da = 0.2, 
     X1a = 0.1,
     X2a = 0.02;

class PField: public FDField {
public: 
    PField(string n, int d[3], double deltas[3]): FDField(n,d,deltas) { }     
    double reactions(std::unordered_map<string,double> params, std::unordered_map<string,FDField*> fields, int x, int y, int z){
        double p = (*values)(x,y,z);
        double a = (*fields["A"])(x,y,z);
        double n = (*fields["N"])(x,y,z);
        return betap*a*n - mp*p; 
    }
    double diffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){ 
        return classicDiffusion(params, fields, Dp, x, y, z);
    }
};

class AField: public FDField {
public:
    AField(string n, int d[3], double deltas[3]): FDField(n,d,deltas) { }     
    double reactions(std::unordered_map<string,double> params, std::unordered_map<string,FDField*> fields, int x, int y, int z){
        return 0; 
    }
    double diffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){ 
        return classicDiffusion(params, fields, Da, x, y, z);
    }
    double chemotaxis(std::unordered_map<string,double> params, std::unordered_map<string,FDField*> fields, int x, int y, int z){
        return classicChemotaxis(params, fields, "A", "P", X1a, x, y, z) + classicChemotaxis(params, fields, "A", "C", X2a, x, y, z);
    }
};

class NField: public FDField {
public:
    NField(string n, int d[3], double deltas[3]): FDField(n,d,deltas) { }    
};

class CField: public FDField {
public:
    CField(string n, int d[3], double deltas[3]): FDField(n,d,deltas) { }    
};


int main(int argc, char *argv[]){    
    int d[3] = {5,5,0}; //Size of each spatial dimension 
    double deltas[3] = {0.1,0.1,0.1}; //Size of the spatial discretizations 
    PField *p = new PField("P",d,deltas);
    AField *a = new AField("A",d,deltas);    
    NField *n = new NField("N",d,deltas);
    CField *c = new CField("C",d,deltas);
    p->saveXYZ();

    n->setInitialValue(3,3,0,10);
    n->setInitialValue(3,2.9,0,10);
    n->setInitialValue(2.9,3,0,10);
    n->setInitialValue(3.1,3,0,10);
    n->setInitialValue(3,3.1,0,10);
    n->setInitialValue(2.8,3,0,10);
    n->setInitialValue(3,2.8,0,10);

    c->setInitialValue(4,4,0,10);
    c->setInitialValue(4.1,4,0,10);
    c->setInitialValue(4,4.1,0,10);
    c->setInitialValue(4.2,4,0,10);
    c->setInitialValue(4,4.2,0,10);

    a->setInitialValue(4,4,0,100);
    a->setInitialValue(4.1,4,0,100);
    a->setInitialValue(4,4.1,0,100);
    a->setInitialValue(4.2,4,0,100);
    a->setInitialValue(4,4.2,0,100);
    
    PdeSolver *solver = new PdeSolver(50, 0.001);
    solver->saveTime();
    vector<double> times = {0.5,1,2,3,4,5,10,20,30,40,50};
    solver->defineReferenceTimes(times);

    solver->addFDField("P", p);
    solver->addFDField("A", a);
    solver->addFDField("N", n);
    solver->addFDField("C", c);
    solver->solve("plot/"); //Run the simulation       
    delete solver; 
    return 0; 
}