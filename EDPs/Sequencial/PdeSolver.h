#ifndef H_PDESOLVER
#define H_PDESOLVER

#include <cstdlib>
#include <exception>
#include <math.h>
#include <omp.h>
#include <iostream> 
#include <sstream>
#include <string>
#include <tuple> 
#include <vector> 
#include <unordered_map> 
#include <fstream> 

using namespace std; 

#include "Array.h"
#include "MultiArray.h"
#include "FDField.h"

class PdeSolver {        
    int tfinal;
    double dt; 
    std::unordered_map<string,double> params;  
    std::unordered_map<string,FDField*> fields; 
    std::vector<std::string> fieldnames; 
    std::vector<double> refTimes; 
    double error; 
public:
    PdeSolver(int final, double deltat){ 
        tfinal = final, dt = deltat; 
        error = 1.0E-06;
    }
    ~PdeSolver(){
        for (std::string n: fieldnames) {
            FDField *f = fields[n]; 
            delete f; 
        }
    }
    bool isReferenceTime(vector<double> times, double ct){
        for (double t: times)
            if (abs(ct - t) <= error)
                return true; 
        return false; 
    }
    void defineReferenceTimes(std::vector<double> times){
        refTimes = times; 
    }
    void addFDField(std::string n, FDField *f){
        if (f != nullptr) {            
            fields[n] = f;
            fieldnames.push_back(n);
        }
    }
    void addParam(std::string pname, double v){
        params[pname] = v;
    }
    void setParams(std::unordered_map<string,double> p){
        params = p; 
    }
    void saveInitialCondition(std::string fpath){
        for (string n: fieldnames) {
            FDField *f = fields[n];
            f->save(0, fpath);
        }
    }
    void advanceTime(double t, std::string fpath){        
        for (string n: fieldnames) {
            FDField *f = fields[n];
            f->update(params, fields, dt); 
            if (isReferenceTime(refTimes, t))            
                f->save(t, fpath);
        }
    }    
    void solve(std::string fpath){
        saveInitialCondition(fpath);
        double t;
        for (t = dt; t <= tfinal; t += dt) {
            advanceTime(t, fpath);
        }
        savePlotConfig();
    }
    void saveTime(){
        ofstream fp;
        fp.open("t.csv");
        fp << "t" << endl;
        double t; 
        for (t = 0; t <= tfinal; t += dt)
            fp << t << endl;
        fp << t;
        fp.close();
    }
    vector<string> splitString(const string& str, const string& delim){
        vector<string> tokens;
        size_t prev = 0, pos = 0;
        do {
            pos = str.find(delim, prev);
            if (pos == string::npos) pos = str.length();
            string token = str.substr(prev, pos-prev);
            if (!token.empty()) tokens.push_back(token);
            prev = pos + delim.length();
        } while (pos < str.length() && prev < str.length());
        return tokens;
    }

    void savePlotConfig(){
/*#se o número de data files é igual a 1, fileNameTime pode ser vazio 
#se fileNameTime é vazio, a primeira coluna do arquivo de dados é o tempo 
numberOfDataFiles: 1 (Determinado com base no número de fields)
dataNames: V, Ap, Apm, Thn, The, Tkn, Tke, B, Ps, Pl, Bm, IgM, IgG, C, I (Nomes dos fields)
fileNames: fieldname1_Ndim.csv, ... 
fileNameTime: t.csv
fileNameX: x.csv
fileNameY: y.csv 
fileNameZ: z.csv 

separatePlots: true
multipleValuesPerTime: false
scales: log10, log10, normal, log10
units: copies/ml, cells/ml [10], S/CO, S/CO, pg/ml, cells/ml
xAxisLabel: days
xAxisRange: 0,100
yAxisRange: 0,30
legendNames: Virus, Immature APC, Mature APC
plotStyle: line*/
        std::unordered_map<std::string,std::string> plotstr; 
        plotstr["numberOfDataFiles"] = to_string(fieldnames.size());
        stringstream ss;
        for (string n: fieldnames)
            if (n == fieldnames.back())
                ss << n;
            else 
                ss << n << ", ";
        plotstr["dataNames"] = ss.str();
        ss.str("");

        for (string n: fieldnames)
            if (n == fieldnames.back())
                ss << fields[n]->getFileName();
            else 
                ss << fields[n]->getFileName() << ", ";
        plotstr["fileNames"] = ss.str();
        ss.str("");
        plotstr["separatePlots"] = "true";
        plotstr["xAxisLabel"] = "days";
        plotstr["yAxisLabel"] = "cells/ml";
        ss << "0-" << tfinal; 
        plotstr["xAxisRange"] = ss.str();
        plotstr["plotStyle"] = "scatter3d";
        ofstream plotfile; 
        plotfile.open("plotconfig_pde.txt");
        for (auto m: plotstr){
            plotfile << m.first << ": " << m.second << endl;
        }        
        plotfile.close(); 
    }

    void readSimulationConfig(string filename){
        // fstream fin;
        // fin.open(path+filename, ios::in);
        // if (fin.is_open()){	
        //     vector<string> row;
        //     string line; 
        //     int t, i = getNewIndex();
        //     double value;            
        //     dataNames[i] = popName;
            
        //     while (fin) {            
        //         row.clear();
        //         if (!getline(fin, line)) break;
        //         row = splitString(line,",");
        //         if (row.size() <= 1) continue; 

        //         try{
        //             t = stod(row[0]);                
        //         }
        //         catch(exception &e){
        //             columnNames[i].push_back(row[0]);                    
        //         }
        //         for (vector<string>::iterator it = std::next(row.begin()); it != row.end(); ++it) {
        //             try{
        //                 value = stod(*it);
        //                 data[i].push_back(make_pair(t,value));
        //             }
        //             catch(exception &e){
        //                 columnNames[i].push_back(*it);                    
        //             }
        //         }     
        //     }
        //     fin.close();	
        // }
        // else 
        //     cout << "Error on file opening." << endl;
        // }
    }
}; 

#endif 