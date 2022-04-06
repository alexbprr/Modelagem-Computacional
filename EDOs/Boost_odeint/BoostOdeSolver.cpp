#include <cmath>
#include <ctime>
#include <boost/numeric/odeint.hpp>
#include <boost/algorithm/string.hpp>
#include <ostream>
#include <sstream>
#include <unordered_map> 
#include <math.h>
#include <fstream> 
#include <numeric>

using namespace std;
using namespace boost::algorithm;
using namespace boost::numeric::odeint;
std::unordered_map<std::string,double> params;

void odesystem(const std::vector<double> &u , std::vector<double> &dudt , const double /* t */) {
    double B = u[0];
    double N = u[1];
    double rb = params["rb"];
    double lnb = params["lnb"];    
    double sn = params["sn"];
    double mn = params["mn"];
    dudt[0] = rb*B/(1 + B) - lnb*N*B;
    dudt[1] = sn*B*N - mn*N; 
}

class OdeSolver {
private:
    int tfinal;
    double dt; 
    std::string filename;    
    std::vector<double> u; //Variables's vector 
    std::vector<std::string> varNames; //Vector of variables names 
    std::vector<double> refTimes;   
    runge_kutta_cash_karp54<std::vector<double>> stepper;
     
public: 
    OdeSolver(){}
    OdeSolver(int tf, double deltat){ 
        tfinal = tf, dt = deltat;         
    }
    ~OdeSolver(){}    
    void defineVarNames(std::vector<std::string> names){ varNames = names; }
    void defineInitialCondition(std::vector<double> values){ u = values; }
    void defineReferenceTimes(std::vector<double> times){ refTimes = times; }
    void addParam(std::string pname, double v){ params[pname] = v; }
    void printInitialCondition(){
        std::cout << "Initial condition:" << std::endl;
        for (double v: u)
            std::cout << v << std::endl;
    }
    void printParams(){
        std::cout << "Parameters:" << std::endl;
        for (std::pair<std::string,double> p: params)
            std::cout << p.first << " = " << p.second << std::endl;
    }
    void printRefTimes(){
        std::cout << "Reference times:" << std::endl;
        for (double v: refTimes)
            std::cout << v << std::endl;
    }
    bool isReferenceTime(std::vector<double> times, double ct){
        int j = 0;
        for (double t: times) {
            if (abs(ct - t) <= 1.0E-06) {
                //Remove t from times 
                times.erase(times.begin()+j);
                return true; 
            }
            j++;
        }
        return false; 
    }
    void createFile(std::string name){
        std::stringstream ss;
        std::ofstream fp;
        ss << name << ".csv";
        filename = ss.str();
        fp.open(filename);
        fp << "t "; 
        for (std::string n: varNames) {
            fp << ", " << n;
        }
        fp << endl; 
        fp.close();
    }
    void save(double t, std::vector<double> values){
        std::ofstream fp;    
        std::stringstream ss;     
        fp.open(filename, ios::app);
        fp << t; 
        for (double v: values)
            fp << "," << v;        
        fp << endl; 
        fp.close();    
    }  
    void solve(){        
        double t = 0;
        createFile("oderesults");
        save(0,u);
        auto c_stepper = make_controlled(1.E-08, 1.E-08, stepper);         
        for (t = dt; t <= tfinal; t += dt) {
            c_stepper.stepper().do_step(odesystem, u, t, dt);
            if (isReferenceTime(refTimes,t))
                save(t,u);
        }
        if (isReferenceTime(refTimes,t))
            save(t,u);        
    }    
  /*
tfinal = 50 
dt = 0.01
refTimes = 1,10,20,30,40,50
----------------
Bacteria = 10
Neutrophil = 2
---------------
rb = 0.5
lnb = 0.05
sn = 0.1 
*/
    std::vector<std::string> splitString(const std::string& str, const std::string& delim){
        std::vector<string> tokens;
        size_t prev = 0, pos = 0;
        do {
            pos = str.find(delim, prev);
            if (pos == string::npos) pos = str.length();
            std::string token = str.substr(prev, pos-prev);
            if (!token.empty()) tokens.push_back(token);
            prev = pos + delim.length();
        } while (pos < str.length() && prev < str.length());
        return tokens;
    }
    void readTimeConfigFromTxt(std::string fname){
        std::fstream fin;
        fin.open(fname, ios::in);
        if (fin.is_open()){	
            std::vector<std::string> row;
            std::string line;             
            while (fin) {            
                row.clear();
                if (!getline(fin, line)) break;                
                row = splitString(line,"=");
                if (row.size() <= 1) break; 
                for(unsigned int i = 0; i < row.size(); i++)
                    trim(row.at(i));        
                try{                    
                    std::string str = row.at(0);
                    if (str == "tfinal") {                    
                        tfinal = std::stod(row.at(1));
                    }
                    else if (str == "dt") {                        
                        dt = std::stod(row.at(1));
                    }
                    else if (str == "refTimes"){
                        //TODO: Get reference to tfinal 
                        row = splitString(row.at(1),",");
                        for(unsigned int i = 0; i < row.size(); i++) {
                            trim(row.at(i));
                            refTimes.push_back(std::stod(row.at(i)));
                        }
                    }
                }
                catch(exception &e){
                    std::cout << e.what() << std::endl;
                }            
            }
            fin.close();
        }
    }
    void readInitialConditionFromTxt(std::string fname){
        std::fstream fin;
        fin.open(fname, ios::in);
        if (fin.is_open()){	
            std::vector<std::string> row;
            std::string line;             
            while (fin) {            
                row.clear();
                if (!getline(fin, line)) break;                
                row = splitString(line,"=");
                if (row.size() <= 1) break; 
                for(unsigned int i = 0; i < row.size(); i++)
                    trim(row.at(i));                
                try{
                    double v = std::stod(row.at(1));
                    varNames.push_back(row.at(0));
                    u.push_back(v);
                }
                catch(exception &e){
                    std::cout << e.what() << std::endl;
                }            
            }
            fin.close();
        }
    }
    void readParamsFromTxt(std::string fname){
        std::fstream fin;
        fin.open(fname, ios::in);
        if (fin.is_open()){	
            std::vector<std::string> row;
            std::string line;             
            while (fin) {            
                row.clear();
                if (!getline(fin, line)) break;                
                row = splitString(line,"=");                
                if (row.size() <= 1) break; 
                for(unsigned int i = 0; i < row.size(); i++)
                    trim(row.at(i));                
                try{
                    double v = std::stod(row.at(1));
                    addParam(row.at(0), v);                    
                }
                catch(exception &e){
                    std::cout << e.what() << std::endl;
                }            
            }
            fin.close();
        }
    }    
    void savePlotConfig(){
        std::unordered_map<std::string,std::string> plotstr; 
        plotstr["numberOfDataFiles"] = "1";
        stringstream ss;
        for (string n: varNames) {
            if (n == varNames.back())
                ss << n;
            else 
                ss << n << ", ";
        }
        plotstr["dataNames"] = ss.str();
        ss.str("");
        plotstr["fileNames"] = filename;
        plotstr["separatePlots"] = "true";
        plotstr["xAxisLabel"] = "days";
        plotstr["yAxisLabel"] = "cells/ml";
        ss << "0-" << tfinal; 
        plotstr["xAxisRange"] = ss.str();
        plotstr["plotStyle"] = "line";
        ofstream plotfile; 
        plotfile.open("plotconfig_ode.txt");
        for (auto m: plotstr){
            plotfile << m.first << ": " << m.second << endl;
        }        
        plotfile.close();        
    }
}; 

int main(){
    //OdeSolver *odesolver = new OdeSolver(100,0.01);
    OdeSolver *odesolver = new OdeSolver();
    //std::vector<std::string> vars = {"Bacteria", "Neutrophil"};
    //odesolver->defineVarNames(vars);
    //std::vector<double> initialValues = {30,5};
    //odesolver->defineInitialCondition(initialValues);
    //std::vector<double> times = {1,3,5,7,9,10,30,50,80,100};
    odesolver->readTimeConfigFromTxt("timeconfig.txt");
    odesolver->readInitialConditionFromTxt("initialcondition.txt");
    //odesolver->defineReferenceTimes(times);    
    /*odesolver->addParam("rb", 0.5);    
    odesolver->addParam("lnb", 0.05);
    odesolver->addParam("sn", 0.1);*/
    //Obs: Para mais de uma execucao, adicionar os parametros apenas uma vez 
    odesolver->readParamsFromTxt("parameters.txt");
    odesolver->solve();    
    odesolver->savePlotConfig();
    delete odesolver; 
    return 0;
}
