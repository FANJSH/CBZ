#include <cstdlib>
#include "RMforFR.h"

using namespace std;

RMforFR_PolynomialFunction::RMforFR_PolynomialFunction(int i)
{
  PutNumPara(i);
};

RMforFR_PolynomialFunction::RMforFR_PolynomialFunction(string filenamepol)
{
  ReadFilepol(filenamepol);
};

void RMforFR_PolynomialFunction::PutNumPara(int i)
{
  num_para = i;
};

int RMforFR_PolynomialFunction::GetNumPara()
{
  return num_para;
};
    
void RMforFR_PolynomialFunction::PutPolynomial(vector<int> &order_inp, double coef_inp)
{
  // It would be necessary to check the same function does not exist
  // in the registered data.
        
  if(order_inp.size()!=num_para){
      cout<<"# Error in PolynomialFunction::PutPolynomial.\n";
      exit(0);
  };

  order.push_back(order_inp);
  for(int i=0; i<num_para;i++){
     //cout<<"Test:"<<order[i][count]<<"\n";
     //count++;
  };
  coef.push_back(coef_inp);
};

void RMforFR_PolynomialFunction::PutPolynomial(int *order_inp, double coef_inp)
{
  vector<int> order_inp2(num_para);
  for(int i=0;i<num_para;i++){
    order_inp2[i]=order_inp[i];
  };
  PutPolynomial(order_inp2, coef_inp);
};

double RMforFR_PolynomialFunction::GetValue(vector<double> x) // To get output parameter via a set of regression functions
{
        int sz=x.size();
        if(sz!=num_para){
            cout<<"# Error in PolynomialFuncition::GetValue.\n";
            cout<<"# Size of the inputted vector x is inconsistent with [num_para].\n";
            exit(0);
        };
        
        double *xinp=new double[sz];
        for(int i=0;i<sz;i++){
            xinp[i]=x[i];
        };
	
	double ret=GetValue(xinp);
	delete [] xinp;
        return ret;
};
    
double RMforFR_PolynomialFunction::GetValue(double *x)
{
        double ret=0.;
        int num_poly=order.size(); // The number of regression functions [I]
        for(int i=0;i<num_poly;i++){    // [I]
            double tmp=coef[i];
            for(int j=0;j<num_para;j++){  // [L]
                tmp*=pow(x[j],order[i][j]);
//                if(i == 0){
//                    cout<<coef[i]<<" "<<x[j]<<" "<<order[i][j]<<endl;
//                };
            };
            ret+=tmp;
            //cout<<i<<" tmp:"<<tmp<<endl;
        };
        return ret;
};
    
void RMforFR_PolynomialFunction::ShowSelf()  // To print out a set of regression functions
{
        int num_poly=order.size(); // The number of regression functions
	cout<<"#-----------------------------------------\n";
	cout<<"# [REGRESSION FUNCTION INFORMATION]\n";
	cout<<"#   The number of terms           : "<<num_poly<<"\n";
	cout<<"#   The number of input variables : "<<num_para<<"\n";
	cout<<"#\n#  ";
        for(int i=0;i<num_poly;i++){
            cout<<" ";
            if(coef[i]>0.)cout<<"+";
            cout<<coef[i];
            for(int j=0;j<num_para;j++){
                int ord=order[i][j];
                if(ord>0){
                    cout<<"[x"<<j<<"]^"<<ord;
                };
            };
        };
        cout<<"\n";
	cout<<"#-----------------------------------------\n";
};

void RMforFR_PolynomialFunction::ReadFilepol(string filenamepol)
{
        ifstream fin;
        fin.open(filenamepol.data(),ios::in);
        if(fin.fail()){
            cout<<"# Error in PolynomiralFunction::ReadFile.\n";
            cout<<"# File named ["<<filenamepol<<"] cannot be found.\n";
            exit(0);
        };

        int itmp;
        double tmp;

        fin>>itmp;  // The number of input parameters at line 1
        PutNumPara(itmp);

        int num_reg_func;  // The number of regression functions at line 2
        fin>>num_reg_func;

        for(int i=0;i<num_reg_func;i++){
            fin>>tmp; // coefficient
            //cout<<i<<" "<<tmp<<endl;

            vector<int> order_inp;

            for(int j=0;j<num_para;j++){
                fin>>itmp; // order
                order_inp.push_back(itmp);
            };
            PutPolynomial(order_inp, tmp);
        };
        
        fin.close();
};

// +++ RMforFR_PostProcessor
RMforFR_PostProcessor::RMforFR_PostProcessor(int i)
{
  PutNumVib(i);
};

RMforFR_PostProcessor::RMforFR_PostProcessor(string filenamepost)
{
  ReadFilepost(filenamepost);
};

void RMforFR_PostProcessor::PutNumVib(int i) //6
{
        num_vib = i;
};
    
void RMforFR_PostProcessor::Putoutputdata(int i) // 16
{
        output_data = i;
};

void RMforFR_PostProcessor::ReadFilepost(string filenamepost)
{
      //test(num_para,num_vi);

        ifstream fin;
        fin.open(filenamepost.data(),ios::in);
        if(fin.fail()){
            cout<<"# Error in PostProcessor::ReadFile.\n";
            cout<<"# File named ["<<filenamepost<<"] cannot be found.\n";
            exit(0);
        };

        int itmp;
        double tmp;
        string str;

        fin>>itmp;
        Putoutputdata(itmp);
         
         fin>>itmp;
         PutNumVib(itmp);
    
        for(int i = 0; i<num_vib; i++)
        {
            fin>>str;
            //cout<<str<<endl;
            
            vector<double> U_inp;
            for(int j = 0; j<output_data; j++){
                fin>>tmp; // U
                U_inp. push_back(tmp);
            };
            U.push_back(U_inp);
            
            for(int k = 0;k<output_data;k++){
                //cout<<"Ui = "<<U[i][k]<<endl;
            };
        };
         //cout<<"Number of U = "<<U.size()<<endl;
         //cout<<"Number of Ui = "<<U[0].size()<<endl;
         
         fin>>str;
         //cout<<str<<endl;
         
         for(int j = 0; j<num_vib;j++)//num_vi
         {
             fin>>tmp;
             S.push_back(tmp);
         };
         
         for(int i = 0; i<output_data;i++)
         {
             fin>>tmp;
             F.push_back(tmp);
         }

        fin.close();
};
    
void RMforFR_PostProcessor::ShowSelf()  // To print out the shape of basis vector U
{
  cout<<"#\n";
  cout<<"# Information on orthnormal basis vectors spanning OUTPUT parameter space\n";
  cout<<"#\n";
  cout<<"#   Dimension of the OUTPUT parameter space : "<<num_vib<<"\n";
  cout<<"#\n";

  for(int i = 0; i<num_vib; i++){
    cout<<"# The "<<i<<"-th basis"<<endl;
    for(int j = 0; j<output_data; j++){
      cout<<U[i][j]<<endl;
    };
    cout<<"\n\n";
  };
};
    
vector<double> RMforFR_PostProcessor::PostSLV(vector<double> & vb)
{
        vector<double> yout(output_data,0.);             //a-hat[I]
        //cout<<"a = "<<num_vi<<"\n";
        
        for(int i = 0; i<num_vib; i++){
            for(int j = 0; j<output_data; j++){
                //cout<<"U = "<<U[i][j]<<"\n";
                yout[j] += U[i][j] * S[i] * vb[i];
            };
        };
        for(int j = 0; j<output_data; j++){
            yout[j] = yout[j] * F[j];
        };

      return yout;
};

// +++ RMforFR_PreProcessor

RMforFR_PreProcessor::RMforFR_PreProcessor(int i)
{
  PutNumVia(i);
};

RMforFR_PreProcessor::RMforFR_PreProcessor(string filenamepre)
{
  ReadFilepre(filenamepre);
};

void RMforFR_PreProcessor::PutNumVia(int i)
{
        num_via = i;
};
    
void RMforFR_PreProcessor::Putinputdata(int i)
{
        input_data = i;
};
    
void RMforFR_PreProcessor::ReadFilepre(string filenamepre)
{
        ifstream fin;
        fin.open(filenamepre.data(),ios::in);
        if(fin.fail()){
            cout<<"# Error in PreProcessor::ReadFile.\n";
            cout<<"# File named ["<<filenamepre<<"] cannot be found.\n";
            exit(0);
        };

        int itmp;
        double tmp;
        string str;

        fin>>itmp;
        Putinputdata(itmp);
        
        fin>>itmp;
        PutNumVia(itmp);
    
        for(int i = 0; i<num_via; i++)
        {
            fin>>str;
            //cout<<str<<endl;
            
            vector<double> U_inp;
            for(int j = 0; j<input_data; j++){
                fin>>tmp; // U
                U_inp. push_back(tmp);
            };
            U.push_back(U_inp);
            
            for(int k = 0;k<input_data;k++){
                //cout<<"Ui = "<<U[i][k]<<endl;
            };
        };
        //cout<<"Number of U = "<<U.size()<<endl;
        //cout<<"Number of Ui = "<<U[0].size()<<endl;
        
        fin>>str;
        //cout<<str<<endl;
        
        for(int j = 0; j<num_via;j++)//num_vi
        {
            fin>>tmp;
            S.push_back(tmp);
        };
        
        //cout<<"number of S = "<<S.size()<<endl;

        fin.close();
};
    
void RMforFR_PreProcessor::ShowSelf()  // To print out the shape of basis vector U
{
  cout<<"#\n";
  cout<<"# Information on orthnormal basis vectors spanning INPUT parameter space\n";
  cout<<"#\n";
  cout<<"#   Dimension of the INPUT parameter space : "<<num_via<<"\n";
  cout<<"#\n";

  for(int i = 0; i<num_via; i++){
    cout<<"# The "<<i<<"-th basis"<<endl;
    for(int j = 0; j<input_data; j++){
      cout<<U[i][j]<<endl;
    };
    cout<<"\n\n";
  };
  /*  
        for(int i = 0; i<num_via; i++){
            cout<<i<<endl;
            for(int j = 0; j<input_data; j++){
                cout<<U[i][j]<<endl;
            };
        };
        cout<<"\n";
  */
};
    
vector<double> RMforFR_PreProcessor::preSLV(vector<double> & xinp)
{
        double pra;
        vector<double> va(num_via,0.);
        //cout<<"a = "<<num_vi<<"\n";
        for(int i = 0; i<num_via; i++){
            for(int j = 0; j<input_data; j++){
                //cout<<" "<<i<<" "<<j<<" "<<"U = "<<U[i][j]<<endl;
                //cout<<"xinp["<<j<<"] ="<<xinp[j]<<endl;
                //cout<<"S["<<i<<"] ="<<S[i]<<endl;
                pra += U[i][j] * xinp[j];
            };
            va[i] = pra / S[i];
            pra = 0;
        };
        //cout<<va.size()<<endl;

      return va;
};


// +++ RMforFR_LowOrderModel

RMforFR_LowOrderModel::RMforFR_LowOrderModel(int n1, int n2, int n3, int n4)
{
        input_data = n1;
        output_data = n2;
        num_via = n3;
        num_vib = n4;
        func.resize(num_vib);
};

void RMforFR_LowOrderModel::ReadPreProcessor(string filenamepre)
{
        prep.ReadFilepre(filenamepre);
};
    
void RMforFR_LowOrderModel::ReadRegressionModel(vector<string> filenamepol)
{
        for(int i=0;i<num_vib;i++)
        {
            func[i].ReadFilepol(filenamepol[i]);
            int tmp=func[i].GetNumPara();
            //cout<<tmp<<endl;
            //cout<<num_via<<endl;
            if(tmp!=num_via)
            {
                cout<<"# Error in LowOrderModel::ReadRegressionModel.\n";
                cout<<"# The number of input parameters is inconsistent.\n";
                exit(0);
            };
            //func[i].ShowSelf();
        };
};
    
void RMforFR_LowOrderModel::ReadPostProcessor(string filenamepost)
{
        postp.ReadFilepost(filenamepost);
};
    
vector<double> RMforFR_LowOrderModel::GetRegressionResult(vector<double> &xinp)
{
        vector<double> va;
        if(num_via != input_data){
            va = prep.preSLV(xinp);
            //prep.ShowSelf();
        }else{
            va = xinp;
        };
        //for(int i=0;i<num_via;i++)
        //{
            //cout<<"va["<<i<<"] = "<<va[i]<<endl; //test
        //};

        vector<double> vb(num_vib);
        for(int i=0;i<num_vib;i++)
        {
            vb[i] = func[i].GetValue(va);
            //cout<<"vb["<<i<<"] = "<<vb[i]<<endl; //test
        };
        
        vector<double> yout;
        if(num_vib != output_data){
            yout=postp.PostSLV(vb);
        }else{
            yout = vb;
        };
        
        //postp.ShowSelf();

        return yout;
};


// +++ RMforFR_Coolingperiod

RMforFR_Coolingperiod::RMforFR_Coolingperiod(int n1, double n2)
{
        output_data = n1;
        period = n2;
};
    
void RMforFR_Coolingperiod::Readfilecp(string filenamecp)
{
        ifstream fin;
        fin.open(filenamecp.data(),ios::in);
        if(fin.fail()){
            cout<<"# Error in Coolingperiod::ReadFile.\n";
            cout<<"# File named ["<<filenamecp<<"] cannot be found.\n";
            exit(0);
        };
        double tmp;
        int itmp;
        for(int i=0;i<output_data;i++){
            fin>>tmp;
            half.push_back(tmp);
            fin>>itmp;
            daughter.push_back(itmp);
        };
};
    
vector<double> RMforFR_Coolingperiod::Coolingcal(vector<double> coolb)
{
        vector<double> coola(output_data);
        double sum = 0;

#if 0	
        double a = -1. * log(2) * (period-1);

        for(int i = 0; i < output_data; i++){
            if(half[i] != 0.)
            {
                coola[i] = coolb[i] * exp(a / half[i]);
                if(daughter[i]!=-1)coola[daughter[i]] += coolb[i] - coola[i];
            }else{
                coola[i] += coolb[i];
            };
            if(coola[i] < 0)
            {
                coola[i] = 0;
            };
            sum += coola[i];
        };
#endif

	real hl_p1=14.29; // year
	real hl_a1=432.6; // year

	real hl_c2=162.94/365.; // year
	real hl_p8=87.7; // year

	real hl_c3=29.1; // year
	real hl_p9=2.411e4; // year

	real hl_c4=18.11; // year
	real hl_p0=6561.; // year

	real lambda_p1=0.693/hl_p1;
	real lambda_a1=0.693/hl_a1;
	real lambda_c2=0.693/hl_c2;
	real lambda_p8=0.693/hl_p8;
	real lambda_c3=0.693/hl_c3;
	real lambda_p9=0.693/hl_p9;
	real lambda_c4=0.693/hl_c4;
	real lambda_p0=0.693/hl_p0;

	real dt=period-1.;

	for(int i=0;i<output_data;i++){
	  coola[i]=coolb[i];
	};

	{
	coola[5]=coolb[5]*exp(-lambda_p1*dt); // Pu-241

	real factor=lambda_p1/(lambda_a1-lambda_p1);
	coola[7]=(coolb[7]-factor*coolb[5])*exp(-lambda_a1*dt)+factor*coolb[5]*exp((lambda_a1-lambda_p1)*dt); // Am-241
	};

	{
	coola[11]=coolb[11]*exp(-lambda_c2*dt); // Cm-242

	real factor=lambda_c2/(lambda_p8-lambda_c2);
	coola[2]=(coolb[2]-factor*coolb[11])*exp(-lambda_p8*dt)+factor*coolb[11]*exp((lambda_p8-lambda_c2)*dt); // Pu-238
	};

	{
	coola[12]=coolb[12]*exp(-lambda_c3*dt); // Cm-243

	real factor=lambda_c3/(lambda_p9-lambda_c3);
	coola[3]=(coolb[3]-factor*coolb[12])*exp(-lambda_p9*dt)+factor*coolb[12]*exp((lambda_p9-lambda_c3)*dt); // Pu-239
	};
	
	{
	coola[13]=coolb[13]*exp(-lambda_c4*dt); // Cm-244

	real factor=lambda_c4/(lambda_p0-lambda_c4);
	coola[4]=(coolb[4]-factor*coolb[13])*exp(-lambda_p0*dt)+factor*coolb[13]*exp((lambda_p0-lambda_c4)*dt); // Pu-240
	};
	
	for(int i=0;i<output_data;i++){
	  sum+=coola[i];
	};

	//cout<<coolb[5]-coola[5]<<" : "<<coolb[7]-coola[7]<<"\n";

        for(int i = 0; i < output_data; i++){
            coola[i] = coola[i] / sum;
            //cout<<"coola["<<i<<"] = "<<coola[i]<<endl;
        };
        return coola;
};


// +++ RMforFR_Showresult
    
RMforFR_Showresult::RMforFR_Showresult(int n1, int n2)
{
        //0:uo2 1:mox 2:core 3:waste 4:TRUCoolingbefore 5:Puenrichment 6:TRUCoolingafter
        if(n1 < 9){
            X = n1;
        }else{
            cout<<"Error:X is not a number from 0 to 4."<<endl;
            exit(0);
        };
    if(n1 != 5){
        output_data = n2;
    }else{
        output_data = n2 * 3;
    }
};
    
void RMforFR_Showresult::resultout(vector<double> result)
{
        int sz=output_data;

        vector<string> name;
        vector<string> unit;
        
        ifstream fin;
        string tmpname="./LIB/nameunit" + IntToString(X);
        fin.open(tmpname.data(),ios::in);	
        if(fin.fail()){
            cout<<"# Error in resultout::ReadFile.\n";
            cout<<"# File named ["<<tmpname<<"] cannot be found.\n";
            exit(0);
        };
        
        string str;
        fin>>str;
        cout<<str<<endl;
        
        for(int i = 0;i<output_data;i++)
        {
            fin>>str;
            name.push_back(str);
            fin>>str;
            unit.push_back(str);
        };

        for(int i=0;i<sz;i++){
            if(X < 2 && i == 0){
                if(result[i] == 0){
                    cout<<i<<" "<<name[i]<<": BWR"<<endl;
                }else{
                    cout<<i<<" "<<name[i]<<": PWR"<<endl;
                };
            }else{
                cout<<i<<" "<<name[i]<<" "<<result[i]<<" "<<unit[i]<<endl;
            }
        };
};

// +++ RMforFR_allcalculation

void RMforFR_allcalculation::ReadInputDataFromFile(string filenamebb)
{
        ifstream fin;
        fin.open(filenamebb.data(),ios::in);
        if(fin.fail()){
            cout<<"# Error in PreProcessor::ReadFile.\n";
            cout<<"# File named ["<<filenamebb<<"] cannot be found.\n";
            exit(0);
        };

        int itmp;
        double tmp;
        
        fin>>itmp;
        X = itmp;

	bool uo2=true;
	if(X == 0){
	}else if(X == 1){
	  uo2=false;
	}else{
          cout<<"Error:enter a number from 0 to 1."<<endl;
          exit(0);
	};

        fin>>tmp;
        real cooling_year = tmp;

	vector<double> tmp2;
	int maxsize=4;
	if(!uo2)maxsize=6;
        for(int i = 0;i<maxsize;i++)
        {
            double tmp;
            fin>>tmp;
            tmp2.push_back(tmp);
        };

	PutInputData(uo2, cooling_year, tmp2);

#if 0 	
        // (tentative input parameters)
        if(X == 0){ //uo2
            input_data = 4;
            output_data = 16;
            num_via = input_data;
            num_vib = 8;
        }else if(X == 1){ //mox
            input_data = 8;
            output_data = 16;
            num_via = input_data;
            num_vib = 8;
        }else{
            cout<<"Error:enter a number from 0 to 1."<<endl;
            exit(0);
        };

        fin>>tmp;
        Y = tmp;
        
        if(X == 0){
            cout<<"uo2-fuel / Cooling period:"<<Y<<" years"<<endl;
        }else if(X == 1){
            cout<<"mox-fuel / Cooling period:"<<Y<<" years"<<endl;
        };
        
        for(int i = 0;i<input_data;i++)
        {
            double tmp;
            fin>>tmp;
            xinp.push_back(tmp);
        };
#endif
	
};

void RMforFR_allcalculation::PutInputData(bool uo2, double cooling_year, vector<double> &inp)
{
  string fuel_name="uo2";  
  if(uo2){ //uo2
    X=0;
    input_data = 4;
    output_data = 16;
    num_via = input_data;
    num_vib = 7;
  }else{ //mox
    X=1;
    input_data = 8;
    output_data = 16;
    num_via = input_data;
    num_vib = 7;
    fuel_name="mox";
  };

  Y = cooling_year; // cooling years
  //cout<<fuel_name<<"-fuel / Cooling period:"<<Y<<"years"<<endl;

  if(inp.size()!=input_data){
    cout<<"# Error in RMforFR_allcalculation::PutInputData.\n";
    cout<<"# The number of input parameters is inconsistent.\n";
    exit(0);
  };

  xinp=inp;
};
    
vector<double> RMforFR_allcalculation::RegressionModelEstimation(bool print_option)
{
        youtal.clear();
        if(print_option){
            RMforFR_Showresult Xin(X, input_data);
            Xin.resultout(xinp);
        };

        RMforFR_LowOrderModel WR(input_data, output_data, input_data, num_vib);
        vector<string> filename_rm;
            
        for(int i=0;i<num_vib;i++){
            filename_rm.push_back("./LIB/polynomial" + IntToString(X) + "_func" + IntToString(i));
        };
        WR.ReadRegressionModel(filename_rm);
        WR.ReadPostProcessor("./LIB/SVD" + IntToString(X) + "_func");

        coolb=WR.GetRegressionResult(xinp);

        if(print_option){
            RMforFR_Showresult TRUb(4, output_data);
            TRUb.resultout(coolb);
        };
        
        int input_fuel = input_data;
        int output_fuel = 1;
        int fuel_via = input_data;
        int fuel_vib = output_fuel;
    
        RMforFR_LowOrderModel Inout_fuel(input_fuel, output_fuel, fuel_via, fuel_vib); // The number of calculation, before SLV

        vector<string> filename_rmfuel;
        for(int i=0;i<fuel_vib;i++)
        {
            filename_rmfuel.push_back("./LIB/polynomial"+IntToString(7+X)+"_func"+IntToString(i));
        };
        
        Inout_fuel.ReadRegressionModel(filename_rmfuel);
        
        vector<double> youtfuel = Inout_fuel.GetRegressionResult(xinp);

        for(int i = 0; i < youtfuel.size(); i++){
            youtal.push_back(youtfuel[i]);
        };
        
        if(print_option){
            RMforFR_Showresult inout_fuel(X+7, output_fuel);
            inout_fuel.resultout(youtfuel);
        };
    
        RMforFR_Coolingperiod cool(output_data, Y);
        cool.Readfilecp("./LIB/Coolingperiod");
        coola = cool.Coolingcal(coolb);
        
        if(print_option){
            RMforFR_Showresult TRUa(6, output_data);
            TRUa.resultout(coola);
	};
    
        //Pu-enrichment calculation
        int input_pu = 16;
        int output_pu = 1;
        int pu_via = 9;
        int pu_vib = output_pu;
        
        RMforFR_LowOrderModel Puen(input_pu, output_pu, pu_via, pu_vib); // The number of calculation, before SLV

        Puen.ReadPreProcessor("./LIB/SVDin_func");
        vector<string> filename_rmpu;
        for(int i=0;i<pu_vib;i++)
        {
            filename_rmpu.push_back("./LIB/polynomial5_func" + IntToString(i));
        };
        
        Puen.ReadRegressionModel(filename_rmpu);
        
        vector<double> youtpu = Puen.GetRegressionResult(coola);

        for(int i = 0; i < youtpu.size(); i++){
            youtal.push_back(youtpu[i]);
        };
        
        if(print_option){
            RMforFR_Showresult Puenrich(5, output_pu);
            Puenrich.resultout(youtpu);
        };
    
        //Core calculation
        int input_c = 16;
        int output_c = 4;
        int c_via = pu_via;
        int c_vib = 3;
        
        RMforFR_LowOrderModel CORE(input_c, output_c, c_via, c_vib); // The number of calculation

        CORE.ReadPreProcessor("./LIB/SVDin_func");
        vector<string> filename_rmc;
        for(int i=0;i<c_vib;i++)
        {
          filename_rmc.push_back("./LIB/polynomial2_func" + IntToString(i));
        };
        
        CORE.ReadRegressionModel(filename_rmc);
        CORE.ReadPostProcessor("./LIB/SVD2_func");
        
        vector<double> youtc = CORE.GetRegressionResult(coola);

        for(int i = 0; i < youtc.size(); i++){
            youtal.push_back(youtc[i]);
        };
        
        if(print_option){
            RMforFR_Showresult core(2, output_c);
            core.resultout(youtc);
        };
    
        //Waste calculation
        int input_w = 16;
        int output_w = 21;
        int w_via = pu_via;
        int w_vib = 17;
        
        RMforFR_LowOrderModel WASTE(input_w, output_w, w_via, w_vib); // The number of calculation

        WASTE.ReadPreProcessor("./LIB/SVDin_func");
        
        vector<string> filename_rmw;
        for(int i=0;i<w_vib;i++){
          filename_rmw.push_back("./LIB/polynomial3_func" + IntToString(i));
        };
        WASTE.ReadRegressionModel(filename_rmw);
        WASTE.ReadPostProcessor("./LIB/SVD3_func");
        
        vector<double> youtw = WASTE.GetRegressionResult(coola);
        
        for(int i = 0; i < youtw.size(); i++){
            youtal.push_back(youtw[i]);
        };
        if(print_option){
            RMforFR_Showresult waste(3, output_w);
            waste.resultout(youtw);
        };

	// ... TRU ratio in the FR loaded fuel is added
	for(int i = 0; i < coola.size(); i++){
	  youtal.push_back(coola[i]);
	};
	
        return youtal;
};


