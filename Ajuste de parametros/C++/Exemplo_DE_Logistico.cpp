#include "de/DifferentialEvolution.h"

#include <cmath>
#include <ctime>
#include <boost/numeric/odeint.hpp>
#include <map> 
#include <math.h>
#include <cstdlib> 
#include <fstream> 
#include <numeric>

using namespace std;
using namespace de;
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;
state_type u;
runge_kutta_cash_karp54<state_type> stepper; 
const double dt = 0.01;
const double tini = 0;
const double tfinal = 50; 
const int pop_size = 1;
const double tol = 1E-06; 

std::vector<double> reference_times; 
std::vector<double> data;

//Initial condition
double N0;
//Parameters 
double r, K;

void odesystem( const state_type &u , state_type &dudt , const double /* t */ ) {
    double N = u[0];        
    dudt[0] = r*N*(1 - N/K);
}

vector<double> advance(double t, std::vector<double> u){
    stepper.do_step(odesystem, u, t, dt);    
    return u;
}

bool isReferenceTime(vector<double> times, double ct){
    for (auto t: times)
        if (abs(ct - t) <= tol)
            return true; 
    return false; 
}

class TestDE : public de::IOptimizable{
public:
    //Solve system and calculate error  
    double EvaluteCost(std::vector<double> inputs) const override{
        double t, error = 0, sumexact = 0; 
        int i = 0;
        u.reserve(1);
        u.resize(1);  
        u[0] = inputs[0];
        r = inputs[1]; 
        K = inputs[2];
        error = 0;
        sumexact = 0;
        for (t = tini; t <= tfinal; t += dt){
            u = advance(t,u);            
            if (isReferenceTime(reference_times,t)){
                error += (u[0] - data[i])*(u[0] - data[i]); //Soma dos quadrados dos erros
                sumexact += u[0]*u[0]; 
                i++; 
            }
        }
        error = sqrt(error/sumexact); //Erro norma 2     
        return error;
    }

    unsigned int NumberOfParameters() const override{
        return 3;
    }

    //Define bounds for each input 
    std::vector<Constraints> GetConstraints() const override{
        std::vector<Constraints> constr(NumberOfParameters());
        constr[0] = Constraints(1, 100, true);
        constr[1] = Constraints(0.01, 1, true);
        constr[2] = Constraints(1, 200, true);
        return constr;
    }

    static bool terminationCondition(const DifferentialEvolution& de){
        if (de.GetBestCost() <= 0.001)
            return true;
        return false; 
    }
};

int main(){
/*
    numberOfDataFiles = 2; 
    data.reserve(numberOfDataFiles);
    columnNames.reserve(numberOfDataFiles);    
    string popName = "t";
    string filename = "n.csv";
    readCsv(filename,popName);
    printData(0);
    
    popName = "N";
    readCsv(filename,popName);
    printData(1);
    
*/
    fstream tfile, datafile;
    tfile.open("data/t.dat", ios::in);
    datafile.open("data/data.dat", ios::in);
    
    //Para cada populacao, ler os valores dos dados experimentais 
    string line; 
    while (tfile >> line) {
        reference_times.push_back(std::stod(line));
    }
    tfile.close();
    line = "";
    //data is global 
    while (datafile >> line) {
        data.push_back(std::stod(line));
    }
    datafile.close();

    cout << "Tempos " << endl;
    for (double t: reference_times)
        cout << t << endl;

    cout << "Dados " << endl;
    for (double t: data)
        cout << to_string(t) << endl;

    u.reserve(1);
    u.resize(1);
    
    TestDE deInstance; 
    de::DifferentialEvolution de(deInstance, 100, std::time(nullptr), 
    true,nullptr,TestDE::terminationCondition);
    
    de.Optimize(50, true);
    cout << "Melhor individuo: "; 
    for (auto v: de.GetBestAgent())
        cout << v << endl;       
    return 0;
}
