#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <random>

#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFitResultPtr.h"

using namespace std;

ifstream fin ( "cf252data.in" );

int N , Z , A;
double yield;
TH2D *PlotNZY = new TH2D ( "Y = Y(N,Z)" , "HistData" , 1000. , 0. , 150. , 1000. , 0. , 100. );
TH2D *PlotAZY = new TH2D ( "Y = Y(A,Z)" , "HistData" , 140 , 50 , 190. , 55 , 20. , 75. );
TH1D *PlotAY = new TH1D ( "Y = Y(A)" , "HistData" , 100. , 70. , 180. );
TH1D *PlotZY = new TH1D ( "Y = Y(Z)" , "HistData" , 100. , 20. , 80. );
void ReadAndPlotTheDatas()
{
	while ( fin )	
	{
		fin >> N >> Z >> yield;
		A = N + Z;
		PlotNZY -> Fill ( N , Z , yield );
		PlotAZY -> Fill ( A , Z , yield );
		PlotAY -> Fill ( A , yield );
		PlotZY -> Fill ( Z , yield );
	}
	PlotZY -> Draw ();
}
double MyFunction ( double *val , double *par )
{
	double x = val[0];
	double y = val[1];
	double meanZ = 14.+par[11]*x+par[12]*x*x;//+par[13]*x*x*x;
	double Gaus1 = par[0]*exp(-0.5*(x-par[1])*(x-par[1])/par[2]/par[2]);
	double Gaus2 = par[3]*exp(-0.5*(x-par[1]-par[4])*(x-par[1]-par[4])/par[5]/par[5]);
	double Gaus3 = par[3]*exp(-0.5*(x-par[1]+par[4])*(x-par[1]+par[4])/par[5]/par[5]);
	double Gaus4 = par[6]*exp(-0.5*(x-par[1]-par[7])*(x-par[1]-par[7])/par[8]/par[8]);
	double Gaus5 = par[6]*exp(-0.5*(x-par[1]+par[7])*(x-par[1]+par[7])/par[8]/par[8]);
	double ZGaus = par[9]*exp(-0.5*(y-meanZ)*(y-meanZ)/par[10]/par[10]);
return (Gaus1 + Gaus2 + Gaus3 + Gaus4 + Gaus5)*ZGaus;
}
double MyFunctionA ( double *valA , double *parA )
{
	double x = valA[0];
	double Gaus1 = parA[0]*exp(-0.5*(x-parA[1])*(x-parA[1])/parA[2]/parA[2]);
	double Gaus2 = parA[3]*exp(-0.5*(x-parA[1]-parA[4])*(x-parA[1]-parA[4])/parA[5]/parA[5]);
	double Gaus3 = parA[3]*exp(-0.5*(x-parA[1]+parA[4])*(x-parA[1]+parA[4])/parA[5]/parA[5]);
	double Gaus4 = parA[6]*exp(-0.5*(x-parA[1]-parA[7])*(x-parA[1]-parA[7])/parA[8]/parA[8]);
	double Gaus5 = parA[6]*exp(-0.5*(x-parA[1]+parA[7])*(x-parA[1]+parA[7])/parA[8]/parA[8]);
	
return Gaus1 + Gaus2 + Gaus3 + Gaus4 + Gaus5;
}
double MyFunctionZ ( double *valZ , double *parZ )
{
	double x = valZ[0];
	double ZGaus = parZ[0]*exp(-0.5*(x-parZ[1])*(x-parZ[1])/parZ[2]/parZ[2]);

return ZGaus;
}
void cfFissionComparison()
{
	ReadAndPlotTheDatas();
	TF2 *FitFunction = new TF2 ( "FitFunction" , MyFunction , 50. , 190., 36. , 70. , 13 );
	TF1 *FitFunctionA = new TF1 ( "FitFunctionA" , MyFunctionA , 50. , 190. , 9 );
	TF1 *FitFunctionZ = new TF1 ( "FitFunctionZ" , MyFunctionZ , 36. , 48.5 , 3 );
	
	FitFunctionA -> SetParameter ( 0 , 1.56908e-05 );
	FitFunctionA -> SetParameter ( 1 , 1.24118e+02 );
	FitFunctionA -> SetParameter ( 2 , 1.61444e+01 );
	FitFunctionA -> SetParameter ( 3 , 6.68403e-02 );
	FitFunctionA -> SetParameter ( 4 , 1.02274e+01 );
	FitFunctionA -> SetParameter ( 5 , 4.45116e+00 );
	FitFunctionA -> SetParameter ( 6 , 1.30728e-02 );
	FitFunctionA -> SetParameter ( 7 , 2.65947e+01 );
	FitFunctionA -> SetParameter ( 8 , 6.22768e+00 );
	PlotAY -> Fit ( FitFunctionA , "REM" );
	
	FitFunctionZ -> SetParameter ( 0 , 1.29031e-01 );
	FitFunctionZ -> SetParameter ( 1 , 4.24164e+01 );
	FitFunctionZ -> SetParameter ( 2 , 2.55416e+00 );
	//PlotZY -> Fit ( FitFunctionZ , "REM" );
	
	FitFunction -> SetParameter ( 0 , 1.56908e-05 );
	FitFunction -> SetParameter ( 1 , 1.21118e+02 );
	FitFunction -> SetParameter ( 2 , 1.61444e+01 );
	FitFunction -> SetParameter ( 3 , 6.68403e-01 );
	FitFunction -> SetParameter ( 4 , 1.62274e+01 );
	FitFunction -> SetParameter ( 5 , 4.46116e+00 );
	FitFunction -> SetParameter ( 6 , 1.30728e-01 );
	FitFunction -> SetParameter ( 7 , 2.65947e+01 );
	FitFunction -> SetParameter ( 8 , 4.82768e+00 );
	FitFunction -> SetParameter ( 9 , 1.29031e-01 );
	FitFunction -> SetParameter ( 10 , 2.55416e+00 );
	FitFunction -> SetParameter ( 11 , 0.253 );
	FitFunction -> SetParameter ( 12 , 5.20758e-05 );
	//FitFunction -> SetParameter ( 13 , 3.74568e-06 );
	//FitFunction -> Draw ( "COLZ" );
	//TFitResultPtr fitTheAZY = PlotAZY -> Fit ( FitFunction , "REM" );
	/*PlotAZY -> Draw ("surf");
	
	//JUST VISUAL PURPOUSES , THE FIT FUNCTION SHOULD BE APPLIED TO PLOT AZY.s
	TH2D* hist = new TH2D( "hist" , "hist" , 110 , 50. , 190. , 50 , 20., 75. );
  	double rA=0. , rZ=0.;
  	for(int iR=0; iR<10000; iR++)
  	{
    		FitFunction->GetRandom2(rA, rZ);
  		hist->Fill(rA, rZ);
 	}	
  	//hist->Draw("surf2");*/
  	
	/*auto FitResult = PlotAZY -> Fit ( FitFunction , "S" );
	double ChiSquared = FitResult -> Chi2();
	int NDF = FitResult -> Ndf();
	cout << '\n' <<"The Value Of ndf = " << NDF << '\n';
	cout << '\n' <<"The Value Of Chi2 = " << ChiSquared << '\n' << '\n';*/
		
}

