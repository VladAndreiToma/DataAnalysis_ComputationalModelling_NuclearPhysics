#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
using namespace std;
ifstream fin ( "cf252data.in" );

int N , Z , A;
double Yield;
TH2D *GeneralHisto = new TH2D ( "histo data" , "General Histo" , 140 , 50. , 190. , 60 , 20. , 80. );

void ReadingData()
{
    while ( fin )
    {
        fin >> N >> Z >> Yield;
        A  = N + Z;
        GeneralHisto -> Fill ( A , Z , Yield );
    }
    GeneralHisto -> Draw ("surf2");
}
void GetTheSlopes()
{
	// note that 1 means first gaussian and 2 the second
	int option;
	cout << "Your option: ";
	cin >> option;
	/*if ( option == 1 )
	{
		TF1 *f1 = new TF1 ( "f1" , "[0]*x+[1]" , 50. , 125. );
		f1 -> SetParameter ( 0 , 0.419994 );
		f1 -> SetParameter ( 1 , 0. );
		GeneralHisto -> Fit ( f1 , "REM" );
	}
	if ( option == 2)
	{
		TF1 *f2 = new TF1 ( "f2" , "[0]*x+[1]" , 125. , 180. );
		f2 -> SetParameter ( 0 , 0.404724 );
		f2 -> SetParameter ( 1 , 0. );
		GeneralHisto -> Fit ( f2 , "REM" );
	}*/
	//TF1 *totalFit = new TF1 ( "totalFit" , "[0]*x+[1]" , 50. , 190. );
	//GeneralHisto -> Fit ( totalFit , "REM" );
}
void improving()
{
    ReadingData();
    GetTheSlopes();
}
