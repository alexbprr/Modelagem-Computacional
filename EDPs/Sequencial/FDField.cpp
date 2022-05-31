#include "FDField.h"
#include "IntMultiArray.h"
#include "MultiArray.h"

FDField::FDField(string n, int d[3], double deltas[3]) {        
    try {
        name = n;
        Ndim = d3D;            
        if (d[1] == 0 && d[2] == 0)
            Ndim = d1D;
        else if (d[2] == 0)
            Ndim = d2D;

        dim[0] = d[0], dim[1] = d[1], dim[2] = d[2];
        dx = deltas[0], dy = deltas[1], dz = deltas[2]; 
        Nx = (d[0]/deltas[0]) + 1, Ny = (d[1]/deltas[1]) + 1, Nz = (d[2]/deltas[2]) + 1;
        std::cout << "Nx " << Nx << std::endl;
        std::cout << "Ny " << Ny << std::endl;
        std::cout << "Nz " << Nz << std::endl;

        xarray = new Array(Nx);
        yarray = new Array(Ny);
        zarray = new Array(Nz);

        int cont = 0;
        for (double x = 0; x <= dim[0]; x += dx) {
            xarray->setValue(cont,x);
            cont++;
        }
        cont = 0; 
        for (double y = 0; y <= dim[1]; y += dy) {
            yarray->setValue(cont,y);
            cont++;
        }
        cont = 0;
        for (double z = 0; z <= dim[2]; z += dz) {
            zarray->setValue(cont,z);
            cont++;
        }        
        int npts[3] = {Nx,Ny,Nz};
        values = new MultiArray(npts);
        nextvalues = new MultiArray(npts);
        sources = new IntMultiArray(npts);
        sinks = new IntMultiArray(npts);
    }
    catch(exception &e){
        cout << e.what() << endl; 
    }
}

FDField::~FDField(){
    delete xarray;
    delete yarray;
    delete zarray;     
    delete values;
    delete nextvalues; 
    delete sources;
    delete sinks; 
}

bool FDField::checkErrors(double v){
    if (isinf(v) || isnan(v)) 
        return true;
    else 
        return false;     
}
bool FDField::convertPosToIndex(double x, double y, double z, int *index){
    for (int i = 0; i < Nx; i++) {
        if (abs ((*xarray)(i) - x) < 10E-08) {
            for (int j = 0; j < Ny; j++) {
                if (abs ((*yarray)(j) - y) < 10E-08)
                    for (int k = 0; k < Nz; k++)
                        if (abs ((*zarray)(k) - z) < 10E-08) {                                                     
                            index[0] = i, index[1] = j, index[2] = k;                            
                            return true;
                        }                      
            }
        }
    }
    return false;    
}
double FDField::getValue(double x, double y, double z){
    int *index = new int[3]; 
    if (convertPosToIndex(x,y,z,index)) {
        return (*values)(index[0],index[1],index[2]); 
    }
    return -1;
}
void FDField::setInitialValue(double x, double y, double z, double v){
    int *index = new int[3]; 
    if (convertPosToIndex(x,y,z,index)) {
        (*values)(index[0],index[1],index[2]) = v; 
    }
}
void FDField::addSource(double x, double y, double z){
    int *index = new int[3]; 
    if (convertPosToIndex(x,y,z,index)) 
        (*sources)(index[0], index[1], index[2]) = 1;
}
int FDField::source(double x, double y, double z){  
    int *index = new int[3]; 
    if (convertPosToIndex(x,y,z,index)) 
        return (*sources)(index[0], index[1], index[2]);
    return 0; 
}
int FDField::hasSource(int x, int y, int z){  
    return (*sources)(x,y,z);    
}
void FDField::addSink(double x, double y, double z){
    int *index = new int[3]; 
    if (convertPosToIndex(x,y,z,index)) 
        (*sinks)(index[0], index[1], index[2]) = 1;       
}
int FDField::sink(double x, double y, double z){  
    int *index = new int[3]; 
    if (convertPosToIndex(x,y,z,index)) 
        return (*sinks)(index[0], index[1], index[2]);
    return 0; 
}
void FDField::printSources(){
    cout << "Source points: " << endl;
    for (int z = 0; z < Nz; z++) 
        for (int y = 0; y < Ny; y++) 
            for (int x = 0; x < Nx; x++) 
                cout << (*xarray)(x) << "," << (*yarray)(y) << "," << (*zarray)(z) << ":" << (*sources)(x,y,z) << endl;
}
void FDField::printSinks(){

}

double FDField::nonlinearDiffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, double D, int x, int y, int z){
    return 0;
}

double FDField::classicDiffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, 
    double D, int x, int y, int z){
    double diffx = 0;
    if (x > 0 && x < Nx-1)
        diffx = (1/(dx*dx))*((*values)(x+1,y,z) - 2*(*values)(x,y,z) + (*values)(x-1,y,z));
    else if (x == 0)
        diffx = (1/(dx*dx))*((*values)(x+1,y,z) - (*values)(x,y,z));            
    else if (x == Nx-1)
        diffx = (1/(dx*dx))*(-(*values)(x,y,z) + (*values)(x-1,y,z));

    double diffy = 0;         
    if (Ndim == d2D || Ndim == d3D) {        
        if (y > 0 && y < Ny-1)
            diffy = (1/(dy*dy))*((*values)(x,y+1,z) - 2*(*values)(x,y,z) + (*values)(x,y-1,z));
        else if (y == 0)
            diffy = (1/(dy*dy))*((*values)(x,y+1,z) - (*values)(x,y,z));            
        else if (y == Ny-1)
            diffy = (1/(dy*dy))*(-(*values)(x,y,z) + (*values)(x,y-1,z));
    }

    double diffz = 0;
    if (Ndim == d3D){
        if (z > 0 && z < Nz-1)
            diffz = (1/(dz*dz))*((*values)(x,y,z+1) - 2*(*values)(x,y,z) + (*values)(x,y,z-1));
        else if (z == 0)
            diffz = (1/(dz*dz))*((*values)(x,y,z+1) - (*values)(x,y,z));            
        else if (z == Nz-1)
            diffz = (1/(dz*dz))*(-(*values)(x,y,z) + (*values)(x,y,z-1));
    }
    return D*(diffx + diffy + diffz);
}

double FDField::diffusion(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){ 
    return 0;
}

double FDField::classicChemotaxis(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, 
    std::string u, std::string attractant, double X, int x, int y, int z){
    double resx = 0, resy = 0, resz = 0, flux_left = 0, flux_right = 0;
    
    FDField *uf = fields[u];
    double n = (*uf->values)(x,y,z);
    double np1x = (*uf->values)(x+1,y,z);
    double nm1x = (*uf->values)(x-1,y,z);
    double np1y = (*uf->values)(x,y+1,z);
    double nm1y = (*uf->values)(x,y-1,z);
    double np1z = (*uf->values)(x,y,z+1);
    double nm1z = (*uf->values)(x,y,z-1);

    FDField *cfield = fields[attractant];
    double c = (*cfield->values)(x,y,z);
    double cp1x = (*cfield->values)(x+1,y,z);
    double cm1x = (*cfield->values)(x-1,y,z);
    double cp1y = (*cfield->values)(x,y+1,z);
    double cm1y = (*cfield->values)(x,y-1,z);
    double cp1z = (*cfield->values)(x,y,z+1);
    double cm1z = (*cfield->values)(x,y,z-1);

    if (x > 0) {
        if(c - cm1x > 0)            
            flux_left = -(c - cm1x)* nm1x / dx;
        else            
            flux_left = -(c - cm1x)* n / dx;
    }
    if(x < Nx-1) {
        if(cp1x - c > 0)        
            flux_right = (cp1x - c)* n / dx;
        else
            flux_right = (cp1x - c)* np1x / dx;
    }    
    resx = (flux_left + flux_right)/dx;

    flux_left = 0, flux_right = 0;
    if (Ndim == d2D || Ndim == d3D) { 
        if (y > 0) {
            if(c - cm1y > 0)            
                flux_left = -(c - cm1y)* nm1y / dy;
            else            
                flux_left = -(c - cm1y)* n / dy;
        }
        if(y < Ny-1) {
            if(cp1y - c > 0)        
                flux_right = (cp1y - c)* n / dy;
            else
                flux_right = (cp1y - c)* np1y / dy;
        }    
        resy = (flux_left + flux_right)/dy;
    }

    flux_left = 0, flux_right = 0;
    if (Ndim == d3D) { 
        if (z > 0) {
            if(c - cm1z > 0)            
                flux_left = -(c - cm1z)* nm1z / dz;
            else            
                flux_left = -(c - cm1z)* n / dz;
        }
        if(z < Nz-1) {
            if(cp1z - c > 0)        
                flux_right = (cp1z - c)* n / dz;
            else
                flux_right = (cp1z - c)* np1z / dz;
        }    
        resz = (flux_left + flux_right)/dz;
    }

    return X*(resx + resy + resz);
}

double FDField::chemotaxis(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){
    return 0; 
}

double FDField::reactions(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z){
    return 0; 
}

bool FDField::update(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, int x, int y, int z, double dt){
    double newvalue = diffusion(params,fields,x,y,z);
    bool error = checkErrors(newvalue);
    if (error) {
        setValue(x,y,z,-1);
        return false;
    }    
    newvalue -= chemotaxis(params,fields,x,y,z);
    error = checkErrors(newvalue);
    if (error) {
        setValue(x,y,z,-1);
        return false;
    } 
    newvalue += reactions(params,fields,x,y,z); 
    error = checkErrors(newvalue);
    if (error) {
        setValue(x,y,z,-1);
        return false;
    }
    else {
        setValue(x,y,z,(*values)(x,y,z) + newvalue*dt);
        return true;
    }
}

void FDField::update(std::unordered_map<std::string,double> params, std::unordered_map<std::string,FDField*> fields, double dt){
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) { 
                if (update(params,fields,x,y,z,dt) == false) {
                    cout << "An error ocurred during the simulation!" << endl;
                    return; 
                }
            }
        }
    }
    MultiArray *tmp = values;
    values = nextvalues;
    nextvalues = tmp; 
}

void FDField::save(double t, std::string fpath){
    std::ofstream fp;    
    std::stringstream ss;
    std::stringstream ss2;
    ss2 << fpath;
    if (dim[2] > 0)
        ss2 << name << "3d_t";
    else if (dim[1] > 0)
        ss2 << name << "2d_t";
    else 
        ss2 << name << "1d_t";
    ss2 << t << ".csv";
    fp.open(ss2.str(), ios::out);    
    if (dim[2] == 0 && dim[1] == 0){
        ss << "x,value"<<endl;
        for (int x = 0; x < Nx; x++)
            ss << (*xarray)(x) << "," << (*values)(x,0,0) << endl;
    }
    else if (dim[2] == 0){
        ss << "x,y,value"<<endl;
        for (int y = 0; y < Ny; y++)
            for (int x = 0; x < Nx; x++)             
                ss << (*xarray)(x) << "," << (*yarray)(y) << "," << (*values)(x,y,0) << endl;
    }
    else {
        ss << "x,y,z,value"<<endl;
        for (int z = 0; z < Nz; z++)
            for (int y = 0; y < Ny; y++)
                for (int x = 0; x < Nx; x++) 
                    ss << (*xarray)(x) << "," << (*yarray)(y) << "," << (*zarray)(z) << "," << (*values)(x,y,z) << endl;
    } 
    fp << ss.str();
    fp.close();
}