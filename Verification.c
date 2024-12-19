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
#include "TSpectrum2.h"

using namespace std;

ifstream fin ( "cf252data.in" );

int N , Z , A;
double YIELDS [1001];
int iterat = 0;
double yield;
//TH2D *PlotNZY = new TH2D ( "Y = Y(N,Z)" , "HistData" , 1000. , 0. , 150. , 1000. , 0. , 100. );
TH2D *PlotAZY = new TH2D ( "Y = Y(A,Z)" , "HistData" , 140 , 50. , 190. , 55 , 20. , 75. );
//TH1D *PlotAY = new TH1D ( "Y = Y(A)" , "HistData" , 100. , 70. , 180. );
//TH1D *PlotZY = new TH1D ( "Y = Y(Z)" , "HistData" , 100. , 20. , 80. );
void ReadAndPlotTheDatas()
{
	while ( fin )	
	{
		fin >> N >> Z >> yield;
		YIELDS [ iterat ] = yield;
		iterat ++;
		A = N + Z;
		//PlotNZY -> Fill ( N , Z , yield );
		PlotAZY -> Fill ( A , Z , yield );
		//PlotAY -> Fill ( A , yield );
		//PlotZY -> Fill ( Z , yield );
	}
}
double MyFunction ( double *val , double *par )
{
	double x = val[0];
	double y = val[1];
	double meanZ = 0 ;
	if( x >= 50 && x < 125 )
		 meanZ = 0.419994*x-1.80108;
	if( x >= 125 && x <= 190 )
		 meanZ = 0.404724*x-1.69983;
	double Gaus0 = par[0]*exp(-0.5*(x-par[1])*(x-par[1])/par[2]/par[2]);
	double Gaus1 = par[3]*exp(-0.5*(x-par[1]-par[4])*(x-par[1]-par[4])/par[5]/par[5]);
	double Gaus2 = par[3]*exp(-0.5*(x-par[1]+par[4])*(x-par[1]+par[4])/par[5]/par[5]);
	double Gaus3 = par[6]*exp(-0.5*(x-par[1]-par[7])*(x-par[1]-par[7])/par[8]/par[8]);
	double Gaus4 = par[6]*exp(-0.5*(x-par[1]+par[7])*(x-par[1]+par[7])/par[8]/par[8]);
	double Gaus5 = par[9]*exp(-0.5*(x-par[1]-par[10])*(x-par[1]-par[10])/par[11]/par[11]);
	double Gaus6 = par[9]*exp(-0.5*(x-par[1]+par[10])*(x-par[1]+par[10])/par[11]/par[11]);
	double ZGaus = 1.*exp(-0.5*(y-meanZ)*(y-meanZ)/par[12]/par[12]);
return (Gaus1 + Gaus2 + Gaus3 + Gaus4 + Gaus5 + Gaus6 + Gaus0)*ZGaus;
}
void Verification()
{
	ReadAndPlotTheDatas();
	TF2 *FitFunction = new TF2 ( "FitFunction" , MyFunction , 50. , 190., 20. , 75. , 14 );
	FitFunction -> SetParameter ( 0 , 1.56908e-06 );
	FitFunction -> SetParameter ( 1 , 1.25018e+02 );
	FitFunction -> SetParameter ( 2 , 1.1444e+00 );
	FitFunction -> SetParameter ( 3 , 4.50403e-02 );
	FitFunction -> SetParameter ( 4 , 1.70004e+01 );
	FitFunction -> SetParameter ( 5 , 4.35116e+00 );
	FitFunction -> SetParameter ( 6 , 1.100728e-02 );
	FitFunction -> SetParameter ( 7 , 2.70947e+01 );
	FitFunction -> SetParameter ( 8 , 3.1568e+00 );
	FitFunction -> SetParameter ( 9 , 1.50031e-02 );
	FitFunction -> SetParameter ( 10 , 9.3000e+00 );
	FitFunction -> SetParameter ( 11 , 2.3668e+00 );
	FitFunction -> SetParameter ( 12 , 0.7000e+00 );
	FitFunction -> SetParameter ( 13 , 1.0000e-04 );
	TH2D* hist = new TH2D( "hist" , "hist" , 140 , 50. , 190. , 55 , 20., 75. );
  	/*double rA=0. , rZ=0.;
  	for(int iR=0; iR<10000; iR++)
  	{
    		FitFunction->GetRandom2(rA, rZ);
  		hist->Fill(rA, rZ);
 	}	
	auto FitResult = PlotAZY -> Fit ( FitFunction , "S" );
	double ChiSquared = FitResult -> Chi2();
	int NDF = FitResult -> Ndf();
	FitResult -> Draw("surf");
	cout << '\n' <<"The Value Of ndf = " << NDF << '\n';
	cout << '\n' <<"The Value Of Chi2 = " << ChiSquared << '\n';
	cout << "The Chi^2 / # of degrees of freedom is: " <<  (double)( ChiSquared / NDF ) << '\n';
	TH2D *ResultHisto = new TH2D ("ratio" , "ratio" , 140 , 50. , 190. , 55 , 20. , 75. );
	TCanvas *c1 = new TCanvas ( "c1" , "c1" , 1920 , 1080 );
	c1 -> Divide (  2 , 2 );
	c1 -> cd(1); 
	PlotAZY -> Draw("COLZ");
	c1 -> cd(2);
	hist -> Draw("surf2");
	ResultHisto = (TH2D*) PlotAZY -> Clone();
	ResultHisto -> GetXaxis() -> SetTitle ("A");
	ResultHisto -> GetYaxis() -> SetTitle ("Z");
	ResultHisto -> SetTitle ("Fission Yield Plot  /  FitFunction Plot");
	ResultHisto -> Divide (hist);
	c1 -> cd(3);
	ResultHisto -> Draw("COLZ");
	c1 -> cd(4);
	PlotAZY -> Draw ("surf2");*/
	PlotAZY -> Draw ("COLZ");
			
}

