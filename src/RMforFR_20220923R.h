#ifndef RMFORFR
#define RMFORFR

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#include "Numeric.h"

using namespace std;

class RMforFR_PolynomialFunction{

  // Input parameter vector x is represented as
  //
  //   x = (x_1, x_2, ..., x_L)
  //
  // Output parameter g(x) is calculated from the following regression model:
  //
  //   g(x) = \sum{i=1}{I} a_i f_i(x)
  //
  // Regression functions f_i(x) is presented as follows:
  //
  //   f_i(x) = \prod{l=1}{L} x_l^{n_{i,l}}
  //

protected:
    int num_para;                 // Number of input parameters [L]
    vector< vector<int> >  order; // List of polynomial order (n_{i,l}): [I][L]
    vector<double> coef;          // Coefficients of regression function (a_i): [I]
    
public:
    RMforFR_PolynomialFunction(int i);
    RMforFR_PolynomialFunction(string filenamepol);
    RMforFR_PolynomialFunction(){};
    ~RMforFR_PolynomialFunction(){};

    void PutNumPara(int i);
    int GetNumPara();
    void PutPolynomial(vector<int> &order_inp, double coef_inp);
    void PutPolynomial(int* order_inp, double coef_inp);    
    double GetValue(vector<double> x); // To get output parameter via a set of regression functions
    double GetValue(double *x); // To get output parameter via a set of regression functions    
    void ShowSelf();  // To print out a set of regression functions
    void ReadFilepol(string filenamepol);
};


class RMforFR_PostProcessor{

protected:
    int num_vib;                   // Number of vi (number of regression models)
    int output_data;                        // number of output parameters
    vector< vector<double> > U;  // [T][num_vi]
    vector<double> S;             // S[I]
    vector<double> V;             // V[I]
    vector<double> F;              // F[I]

public:
    RMforFR_PostProcessor(int i);
    RMforFR_PostProcessor(string filenamepost);
    RMforFR_PostProcessor(){};
    ~RMforFR_PostProcessor(){};

    void PutNumVib(int i);    
    void Putoutputdata(int i);
    void ReadFilepost(string filenamepost);
    void ShowSelf();  // To print out the shape of basis vector U
    vector<double> PostSLV(vector<double> & vb);
};


class RMforFR_PreProcessor{

protected:
    int num_via;                   // Number of vi (number of regression models)
    int input_data;                        // TRU[16]      (number of output parameters)
    vector< vector<double> > U;  // [T][num_vi]
    vector<double> S;             // S[I]
    vector<double> V;             // V[I]
    
public:
    RMforFR_PreProcessor(int i);
    RMforFR_PreProcessor(string filenamepre);
    RMforFR_PreProcessor(){};
    ~RMforFR_PreProcessor(){};

    void PutNumVia(int i);
    void Putinputdata(int i);
    void ReadFilepre(string filenamepre);
    void ShowSelf();  // To print out the shape of basis vector U
    vector<double> preSLV(vector<double> & xinp);
};


class RMforFR_LowOrderModel{

protected:
    int input_data, output_data;
    int num_via, num_vib;
    vector<RMforFR_PolynomialFunction> func;
    RMforFR_PreProcessor prep;
    RMforFR_PostProcessor postp;
public:
    RMforFR_LowOrderModel(int n1, int n2, int n3, int n4);
    RMforFR_LowOrderModel(int n1, int n2, int n3, bool beoraf);
    ~RMforFR_LowOrderModel(){};

    void Initialize1(int n1, int n2, int n3, int n4);
    void Initialize2(int n1, int n2, int n3, bool beoraf);
    
    void ReadPreProcessor(string filenamepre);
    void ReadRegressionModel(vector<string> filenamepol);
    void ReadPostProcessor(string filenamepost);
    vector<double> GetRegressionResult(vector<double> &xinp);
    void ShowPreProcessor(){prep.ShowSelf();};
    void ShowPostProcessor(){postp.ShowSelf();};
    RMforFR_PolynomialFunction GetRegressionFunction(int i){return func[i];};
};


class RMforFR_Coolingperiod{

protected:
    int output_data;
    double period;
    vector<double> half;
    vector<int> daughter;

public:
    RMforFR_Coolingperiod(int n1, double n2);
    ~RMforFR_Coolingperiod(){};
    
    void Readfilecp(string filenamecp);
    vector<double> Coolingcal(vector<double> coolb);
};


class RMforFR_Showresult{
    
protected:
    int output_data;
    int X;

public:
    RMforFR_Showresult(int n1, int n2);
    ~RMforFR_Showresult(){};
    
    void resultout(vector<double> result);
};

class RMforFR_allcalculation{
    
protected:
    int X;
    double Y;
    int input_data, output_data, num_via, num_vib;
    vector<double> xinp;
    vector<double> coolb;
    vector<double> coola;
    vector<double> youtc;
    vector<double> youtw;
    vector<double> youtal;

public:
    RMforFR_allcalculation(){};
    RMforFR_allcalculation(string filenamebb){ReadInputDataFromFile(filenamebb);};    
    ~RMforFR_allcalculation(){};    

    void ReadInputDataFromFile(string filenamebb);
    void PutInputData(bool uo2, double cooling_year, vector<double> &inp);
    vector<double> RegressionModelEstimation(bool print_option=true);
};

#endif
