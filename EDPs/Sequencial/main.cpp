#include "PdeSolver.h"
#include <random>

/*Use (*xarray)(x),(*yarray)(y) and (*xarray)(z) to access the values of the coordinates of a point
*/ 
//Create a field for each population of the system  
class BField: public FDField {
public:
    //Reactions, diffusion and chemotaxis can be defined here 
    BField(string n, int d[3], double deltas[3]): FDField(n,d,deltas) { }     
    double reactions(std::unordered_map<string,double> params, std::unordered_map<string,FDField*> fields, int x, int y, int z){
        double b = (*values)(x,y,z);
        double n = (*fields["neutrophil"])(x,y,z);
        //return params["Rb"]*b - params["Lnb"]*n*b;
        return 0; 
    }
    double diffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){ 
        return classicDiffusion(params, fields, params["Db"], x, y, z);
    }
};

class NField: public FDField {
public:
    NField(string n, int d[3], double deltas[3]): FDField(n,d,deltas) { }     
    double reactions(std::unordered_map<string,double> params, std::unordered_map<string,FDField*> fields, int x, int y, int z){
        double n = (*values)(x,y,z);
        double c = (*fields["proinflammatory-cytokines"])(x,y,z);        
        double sourceN = 0;
        if (hasSource(x,y,z)){
            sourceN = params["NBlood"]*(params["Pnfix"] + params["Pnrel"]*c/(c + params["Kc"]));
        }        
        //return sourceN - params["mn"]*n;
        return 0;
    }
    double diffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){ 
        return classicDiffusion(params, fields, params["Dn"], x, y, z);
    }
    double chemotaxis(std::unordered_map<string,double> params, std::unordered_map<string,FDField*> fields, int x, int y, int z){
        return classicChemotaxis(params, fields, "neutrophil", "proinflammatory-cytokines", params["Xn"], x, y, z);
    }
};

class CField: public FDField {
public:
    CField(string n, int d[3], double deltas[3]): FDField(n,d,deltas) { }     
    double reactions(std::unordered_map<string,double> params, std::unordered_map<string,FDField*> fields, int x, int y, int z){
        double c = (*values)(x,y,z);
        double b = (*fields["bacteria"])(x,y,z);
        double n = (*fields["neutrophil"])(x,y,z);        
        //return params["Betacn"]*b - params["Mc"]*c;
        return 0;
    }
    double diffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){ 
        return classicDiffusion(params, fields, params["Dc"], x, y, z);
    }
};

int main(int argc, char *argv[]){    
    int d[3] = {5,0,0}; //Size of each spatial dimension 
    double deltas[3] = {0.1,0.1,0.1}; //Size of the spatial discretizations 
    BField *bacteria = new BField("bacteria",d,deltas);    
    NField *neutrophil = new NField("neutrophil",d,deltas);
    string cname = "proinflammatory-cytokines";
    CField *c = new CField(cname,d,deltas);
    bacteria->saveXYZ();

    neutrophil->addSource(1.0, 0, 0);
    neutrophil->addSource(3.0, 0, 0);
    neutrophil->addSource(3.0, 5, 0);
    //neutrophil->printSources();
    
    bacteria->setInitialValue(2.5,0,0,10);
    
    PdeSolver *solver = new PdeSolver(50, 0.01);
    solver->saveTime();
    vector<double> times = {0.5,1,2,3,4,5,10,20,30,40,50};
    solver->defineReferenceTimes(times);

    solver->addFDField("bacteria", bacteria);
    solver->addFDField("neutrophil", neutrophil);
    solver->addFDField(cname, c);

    solver->addParam("Rb", 0.2);
    solver->addParam("Lnb", 0.1);
    solver->addParam("Db", 0.02);
    solver->addParam("Pnfix", 0.0001);
    solver->addParam("Pnrel", 0.2);
    solver->addParam("Kc", 10);
    solver->addParam("NBlood", 1E05);
    solver->addParam("Mn", 0.5);
    solver->addParam("Dn", 0.05);
    solver->addParam("Xn", 0.03);
    solver->addParam("Betacn", 0.7);
    solver->addParam("Mc", 0.7);
    solver->addParam("Dc", 0.1);
    solver->solve("plot/"); //Run the simulation       
    delete solver; 
    return 0; 
}