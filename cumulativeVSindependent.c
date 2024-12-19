#include "MyLib.hh"
#include "THStack.h"
#include "TStyle.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "THistPainter.h"
#include "TAttMarker.h"
using namespace std;
ifstream finC( "cumulativeDataInput.txt" );
ifstream finI( "independentDataInput.txt" );
int Ai , Ac , Aimax = 0 , Aimin = 999999 , Acmin = 999999, Acmax = 0 ;
double Yi[107] , Yc[107];
double ERRi[107] , ERRc[107];
string PerERRi[107] , PerERRc[107];

TH1D* IndepHisto = new TH1D ( "IndepHisto" , "Y = Y(A)" , 120 , 60. , 180. );
TH1D* CumulHisto = new TH1D ( "CumulHisto" , "Y = Y(A)" , 120 , 60. , 180. );
 
void ReadData()
{
	// read the independent spectra Y=Y(A)
	int i = 0;
	while ( finI )
	{
		finI >> Ai >> Yi[i] >> ERRi[i] >> PerERRi[i];
		if ( Ai > Aimax )
			Aimax = Ai;
		if ( Ai < Aimin )
			Aimin = Ai;
		i++;
	}
	int j = 0;
	finI . close(); // closing the independent data set
	while ( finC )
	{
		finC >> Ac >> Yc[j] >> ERRc[j] >> PerERRi[j];
		if ( Ac > Acmax )
			Acmax = Ac;
		if ( Ac < Acmin )
			Acmin = Ac;
		j++;
	}
	finC . close(); // closing the cummulative data set
}
void FillTheHistograms()
{
	for ( int i = Aimin ; i <= Aimax ; i++ )
	{
		IndepHisto -> SetBinContent ( (double)i - 60. , (double)log(Yi[i-66]) );
		IndepHisto -> SetBinError ( (double)i - 60. , (double)(ERRi[i-66]) );
	}
	for ( int j = Acmin ; j <= Acmax ; j++ )
	{
		CumulHisto -> SetBinContent ( (double)j - 60. , (double)log(Yc[j-66]) );
		CumulHisto -> SetBinError ( (double)j - 60. , (double)(ERRc[j-66]) );
	}
}
void PlotData()
{
	//gPad -> SetLogy(1);
	IndepHisto -> SetMarkerColor( kRed );
	CumulHisto -> SetMarkerColor( kBlue );
	IndepHisto -> Draw("PMC");
	CumulHisto -> Draw( "same,PMC" );	
}		
void cumulativeVSindependent()
{
	ReadData();
	FillTheHistograms();
	PlotData();
}
