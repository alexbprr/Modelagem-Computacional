#include <iostream> 
#include <fstream> 
#include <omp.h>
#include <sstream>
#include <random>
#include <string> 
#include <unordered_map>
#include <cmath> 
#include <algorithm> 
#include <boost/algorithm/string.hpp>
using namespace std; 
using namespace boost::algorithm;

std::random_device rd; 
std::mt19937 mersenne_engine {rd()}; // Mersenne Twister
std::uniform_real_distribution<double> distribution (0.0,1.0);

class GillespieSolver {
    int tfinal;    
    double dt;     
    std::string filename;     
    std::vector<int> u; 
    std::vector<std::string> varNames; 
    std::vector<double> refTimes; 
    std::unordered_map<std::string,double> params;

public:
    GillespieSolver(){}
    GillespieSolver(int tf, double deltat){ tfinal = tf, dt = deltat; }
    GillespieSolver(int tf, double deltat, std::vector<std::string> names, std::vector<int> u_, 
        std::vector<double> times){
        tfinal = tf, dt = deltat;
        varNames = names;
        refTimes = times;
        u = u_; 
    }
    ~GillespieSolver(){}
    void defineVarNames(std::vector<std::string> names){ varNames = names; }
    void defineReferenceTimes(std::vector<double> times){ refTimes = times; }
    void defineInitialCondition(std::vector<int> u_){ u = u_;}
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
        for (double t: times)
            if (abs(ct - t) <= 1.0E-06)
                return true; 
        return false; 
    }    
    void saveInitialCondition(std::string name){
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
    void save(double t, std::vector<int> values){
        std::ofstream fp;    
        std::stringstream ss;     
        fp.open(filename, ios::app);
        fp << t; 
        for (double v: values)
            fp << ", " << v;        
        fp << endl; 
        fp.close();    
    }     
    double gillespie() {
        double x, tau, st; 
        double H = u[0];
        double P = u[1];

        double term0_H = params["rh"]*H;
        double term1_P = params["a"]*H*P /(1 + params["tetah"]*H) ;
        double term2_P = params["mp"]*P;
        
        st = term0_H + term1_P + term2_P;

        x = distribution(mersenne_engine);
        tau = -std::log(distribution(mersenne_engine))/st;
                
        //Testing which reaction will occur                 
        if (x < (term0_H)/st) {
            H++;
        }
        else if (x < (term0_H + term1_P)/st) {
            P++;
            H--;
        }
        else if (x < (term0_H + term1_P + term2_P)/st) {
            P--;
        }

        //Updating the values  
        u[0] = H;
        u[1] = P; 
        return tau; 
    }
    void advanceTime(){
        double time = 0.0;
        while(time < dt){
            time += gillespie();
        }
    }
    void solve(){
        double t = 0;
        save(0,u);
        for (t = dt; t <= tfinal; t += dt) {
            advanceTime();
            if (isReferenceTime(refTimes,t))
                save(t,u);
        }
        if (isReferenceTime(refTimes,t))
            save(t,u);         
    }
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
    void savePlotConfigFile(int nRuns, std::vector<std::string> filenames){
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
        for (string s: filenames) {
            if (s == filenames.back())
                ss << s;
            else 
                ss << s << ", ";
        }        
        plotstr["fileNames"] = ss.str();
        ss.str("");
        plotstr["separatePlots"] = "true";
        plotstr["plotManyResults"] = "true";
        plotstr["numberOfResults"] = to_string(nRuns);
        plotstr["xAxisLabel"] = "days";
        plotstr["yAxisLabel"] = "cells/ml";
        ss << "0-" << tfinal; 
        plotstr["xAxisRange"] = ss.str();
        plotstr["plotStyle"] = "line";
        ofstream plotfile; 
        plotfile.open("plotconfig_gillespie.txt");
        for (auto m: plotstr){
            plotfile << m.first << ": " << m.second << endl;
        }
        plotfile.close(); 
    }    
};

int main(){
    int nRuns = 5;
    /*
    std::vector<std::string> vars = {"Prey", "Predator"};
    std::vector<int> initialValues = {30,5};
    double tf = 10000; 
    std::vector<double> times = {1,3,5,7,9,10,15,20,25,30,50,60,70,80,90,100,tf/10.,tf};*/    
    std::vector<std::string> filenames;
    #pragma omp parallel for num_threads(nRuns)
    for (int i = 0; i < nRuns; i++){        
        //GillespieSolver *gsolver = new GillespieSolver(tf,0.1,vars,initialValues,times);
        GillespieSolver *gsolver = new GillespieSolver();
        gsolver->readTimeConfigFromTxt("entrada/timeconfig.txt");
        gsolver->readInitialConditionFromTxt("entrada/initialcondition.txt");
        gsolver->readParamsFromTxt("entrada/parameters.txt");
        if (omp_get_thread_num() == 0) {
            gsolver->printParams();
            gsolver->printInitialCondition();
            gsolver->printRefTimes();
        }
        stringstream filename; 
        filename << "gilles_run_" << i;
        filenames.push_back(filename.str()); 
        gsolver->saveInitialCondition(filename.str());
        /*gsolver->addParam("rh",0.2);
        gsolver->addParam("a",0.05);
        gsolver->addParam("mp",0.1);
        gsolver->addParam("tetah",0.15);*/
        gsolver->solve();
        if (omp_get_thread_num() == 0)
            gsolver->savePlotConfigFile(nRuns, filenames);
        delete gsolver; 
    }
    return 0;
}
