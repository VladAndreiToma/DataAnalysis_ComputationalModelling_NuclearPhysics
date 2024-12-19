

// ******************************** including a self-created library that contains all the needed main libraries from root ***********************************************
#include "MyLib.hh"
using namespace std;  // using the standard namespace for this program
// *************************************************************************************************************************************************************************
ifstream fin ( "DataOutput2D.txt" );  // reading the data from the constructed file .txt , this is the result after a script used to interpret the raw data from NNDC (those were matrix aranged so a little work had to be done to adjust them as column list)
// **************************** the parameters of the code ***************************************************************************************************************
double yield;
int A , Z;
bool option = 0;
int RangeAMax = 200 , RangeAMin = 1 , RangeZMax = 100 , RangeZMin = 1;
TH2D *MainHisto = new TH2D ( "datas" , "Histogram from NNDC Experimental Data - Fission Yields" , 140 , 50. , 190. , 55 , 20. , 75. ); // the MAIN histogram of the entire script
// **************************** end of parameters list ******************************************************************************************************************
// **************************** READING THE DATA FROM THE TXT FILE ************************************************************************************************************************
void ReadData()
{
	int i = 0;
	while ( fin )
	{
		fin >> yield >> A >> Z;
		MainHisto -> SetBinContent ( A - 50 , Z - 20 , yield/100 ); // divide the yields with 100 to get them as percentages
	}
}
//**************************** END OF READING THE FILE ******************************************************************************************************************
// ************************** PLOTING THE DATAS , METHOD OF A PLOT IS USER'S CHOICE *************************************************************************************
void PlotTheData()
{
	if( option )
		MainHisto -> Draw( "surf2" );
	else 
		MainHisto -> Draw( "colz" );
}
// **************************** END OF DATA PLOTTING ***********************************************************************************************************************
// HERE I CONSTRUCTED A FIT FUCNTION TO PROVIDE MATHEMATICAL INTERPRETATION FOR THE FISSION SPECTRA OF YIELDS IN TERMS OF PROTON NUMBER AND ATOMIC MASS NUMBER *************
// Everything has to be done via trial and error tehnique, with a previous analytical interpretation to see what is the most appropriate form of the function, before REM method is called to find the actual values of the entire list for funcParam. *funcParam is an array of parameters and *ranges an array of stepts on the x and y axis. both are pointers to themselves and are dynamically modified throughout the execution of the algorithm .........
double CompoundFitFunction( double* ranges , double *funcParam )
{
	double a = ranges[0];
	double z = ranges[1];
	double theta1 = 0;
	double theta2 = 0;
	double zMean = 0;
	if ( a >= 50 && a < 118 )
    		zMean = 0.416788 * a - 2.00559;
    	if ( a >=128 && a <=190 )
    		zMean = 0.405715 * a - 2.42845;
	if( a >= 118 && a < 128 )
		zMean = 47.5451;
	double Gauss1A = funcParam[0] * exp ( - 0.5 * ( a - funcParam[1] ) * ( a - funcParam[1] ) / funcParam[2] / funcParam[2] );
	double Gauss2A = funcParam[3] * exp ( - 0.5 * ( a - funcParam[1] - funcParam[4] ) * ( a - funcParam[1] - funcParam[4] ) / funcParam[5] / funcParam[5] );
	double Gauss3A = funcParam[3] * exp ( - 0.5 * ( a - funcParam[1] + funcParam[4] ) * ( a - funcParam[1] + funcParam[4] ) / funcParam[5] / funcParam[5] );
	double Gauss4A = funcParam[6] * exp ( - 0.5 * ( a - funcParam[1] - funcParam[7] ) * ( a - funcParam[1] - funcParam[7] ) / funcParam[8] / funcParam[8] );
	double Gauss5A = funcParam[6] * exp ( - 0.5 * ( a - funcParam[1] + funcParam[7] ) * ( a - funcParam[1] + funcParam[7] ) / funcParam[8] / funcParam[8] );
	double Gauss6A = funcParam[9] * exp ( - 0.5 * ( a - funcParam[1] - funcParam[10] ) * ( a - funcParam[1] - funcParam[10] ) / funcParam[11] / funcParam[11] );
	double Gauss7A = funcParam[9] * exp ( - 0.5 * ( a - funcParam[1] + funcParam[10] ) * ( a - funcParam[1] + funcParam[10] ) / funcParam[11] / funcParam[11] );
	double ZGauss = 1.0 * exp ( -0.5 * ( z - zMean ) * ( z - zMean ) / funcParam[12] / funcParam[12] );
	
return ( Gauss1A + Gauss2A + Gauss3A + Gauss4A + Gauss5A + Gauss6A + Gauss7A ) * ZGauss;
}
// **************************** END OF SELF CONSTRUCTED FIT FUNCTION *******************************************************************************************************
// **************************** A CLONE FUNCTION OF MY FIT FUNCTION TO PROVIDE YIELDS TO FILL THE CHECKING HISTOGRAM *******************************************************
// the parameters are not allowed to change, they are pre - established and they remain cosntant thorughout the execution of the algo
double MyGeneratedYield( int genA , int genZ )
{
	double zMean = 0;
	if ( genA >= 50 && genA < 120 )
    		zMean = 0.416788 * genA - 2.00559;
    	if ( genA >= 125 && genA <= 190 )
    		zMean = 0.405715 * genA - 2.42845; 
	if( genA >= 120 && genA < 125 )
		zMean = 47.5451;
	double Gauss1A = 1.35e-02 * exp ( - 0.5 * ( genA - 1.23e+02 ) * ( genA - 1.22e+02 ) / 2.46e-01 / 2.46e-01 );
	double Gauss2A = 4.30e+00 * exp ( - 0.5 * ( genA - 1.23e+02 - 1.80e+01 ) * ( genA - 1.23e+02 - 1.80e+01  ) / 4.72e+00 / 4.72e+00 );
	double Gauss3A = 4.30e+00 * exp ( - 0.5 * ( genA - 1.23e+02 + 1.80e+01 ) * ( genA - 1.23e+02 + 1.80e+01  ) / 4.72e+00 / 4.72e+00 );
	double Gauss4A = 1.10e+00 * exp ( - 0.5 * ( genA - 1.23e+02 - 2.81e+01 ) * ( genA - 1.23e+02 - 2.81e+01 ) / 4.67e+00 / 4.67e+00 );
	double Gauss5A = 1.10e+00 * exp ( - 0.5 * ( genA - 1.23e+02 + 2.81e+01 ) * ( genA - 1.23e+02 + 2.81e+01 ) / 4.67e+00 / 4.67e+00 );
	double Gauss6A = 2.36e+00 * exp ( - 0.5 * ( genA - 1.23e+02 - 11.6e+00  ) * ( genA - 1.23e+02 - 10.6e+00 ) / 2.83e+00 / 2.83e+00 );
	double Gauss7A = 2.36e+00 * exp ( - 0.5 * ( genA - 1.23e+02 + 11.6e+00  ) * ( genA - 1.23e+02 + 10.6e+00 ) / 2.83e+00 / 2.83e+00 );
	double ZGauss = 1.0 * exp ( -0.5 * ( genZ - (int)zMean ) * ( genZ - (int)zMean ) / 1.00e+00 / 1.00e+00 );
	
return ( Gauss1A + Gauss2A + Gauss3A + Gauss4A + Gauss5A + Gauss6A + Gauss7A ) * ZGauss;
}
// ******************* THE END OF THE YIELD FILL FUNCTION **********************************************************************************************************************
// ******************* THE MAIN CODE WHERE EVERYTHING IS EXECUTED AND CALLED ***************************************************************************************************
void UsrDefBin2D()
{
	// lines that allow the user to define the option:
	cout << "*******GIVE YOUR OPTION , 0 OR 1 , 0 IS FOR COLZ PLOT , 1 IS FOR SURF2 PLOT ..... OPTION = ";
	cin >> option; 
	// ***************** SOME GUIDING FUNCTION FOR THE DEPENDENCE: ZMEAN = ZMEAN (A) FOR CF SOURCE ******************************************
	TF1* f1 = new TF1 ( "f1" , "[0]+[1]*x" , 50. , 119. );
	f1 -> SetParameter ( 0 , -2. );
	f1 -> SetParameter ( 1 , 0.5 );
	TF2* f2 = new TF2 ( "f2" , "pol1" , 119. , 123. );
	//***************************************************************************************************************************************
	TCanvas* MyCanvas = new TCanvas ( "c1" , "Plots" , 1920. , 1080. );   // canvas object defined and instantiated
	MyCanvas -> Divide ( 2 , 2 );   // division in a 2*2 format
	MyCanvas -> cd ( 1 );    // cell of indices: 1 , 1 
	ReadData();     // reading
	PlotTheData();  // visual interpretations
	MyCanvas -> cd ( 2 );   // cell of indices 1 , 2
	TH2D* CloneMainHisto = (TH2D*) MainHisto -> Clone();   // making a clone of my main histogram for later use
	TF2* MyFitFunction = new TF2 ( "MyFitFunction" , CompoundFitFunction , 50. , 190. , 20. , 75. , 13 );   // defining my fit function whith the compound fit function properties and parameters
	// setting the parameters of the fit function , from here REM finds the exact values , the values that i wrote for the parameters are in the wanted range but not exactly the most precise ones
	MyFitFunction -> SetParameter ( 0 , 1.35e-02 );
	MyFitFunction -> SetParameter ( 1 , 1.23e+02 );
	MyFitFunction -> SetParameter ( 2 , 2.46e-01 );
	MyFitFunction -> SetParameter ( 3 , 4.30e+00 );
	MyFitFunction -> SetParameter ( 4 , 1.80e+01 );
	MyFitFunction -> SetParameter ( 5 , 4.72e+00 );
	MyFitFunction -> SetParameter ( 6 , 1.10e-01 );
	MyFitFunction -> SetParameter ( 8 , 4.97e+00 );
	MyFitFunction -> SetParameter ( 9 , 2.36e+00 );
	MyFitFunction -> SetParameter ( 10 , 11.6e+00 );
	MyFitFunction -> SetParameter ( 7 , 2.81e+01 );
	MyFitFunction -> SetParameter ( 11 , 3.03e+00 );
   	MyFitFunction -> SetParameter ( 12 , 1.00e+00 );
   	// ******************* end of setting the function parameters ***************************888
   	int Runs = 100000;
   	double Zrand = 0 , Arand = 0;
   	TH2D *GeneratedHisto = new TH2D( "Generated Histo" , "Y = Y (A,Z)" , 140 , 50. , 190. , 55 , 20. , 75. );  // generate a histogram from my fit function
   	while( Runs )
   	{
   		MyFitFunction -> GetRandom2 ( Arand , Zrand ); // get random A , Z from the allowed domain of MyFitFunction 
   		// A and Z are sent as addresses , modified and then returned back;
   		GeneratedHisto -> SetBinContent ( (int)Arand - 50 , (int)Zrand - 20 , MyGeneratedYield ( (int)Arand , (int)Zrand )); // now the yield is calculated with the clone function: #64
   		// everything is setted bin by bin equidistant, to respect the scallation of the whole histogram
   		Runs --; // the filling process is done Runs times;
   	}   	
	if ( option )
		GeneratedHisto -> Draw ( "surf2" );
	else
		GeneratedHisto -> Draw( "colz" );
	MyCanvas -> cd ( 3 );   // cell of indices: 2 , 1
	TH2D* RatioOfHistos = new TH2D(*MainHisto); RatioOfHistos -> Divide ( GeneratedHisto ); // here it is checked the resamblance of the 2 histos by division
	// best values need to bee around 1 or even 1 for perfection, that means that the fit reconstructs very well the datas and provide outcomes where raw experimental data are not measured (because it is not statistically that possible, or from other reasons ).
	RatioOfHistos -> SetNameTitle ( "RatioOfHistos" , "DataHisto/FitFunctionHisto");
	if ( option )
		RatioOfHistos -> Draw ( "surf2" );
	else
		RatioOfHistos -> Draw ( "colz" );
	MyCanvas -> cd ( 4 );   // the cell of indices: 2 , 2
	TFitResultPtr FitResult = CloneMainHisto -> Fit ( MyFitFunction , "REMS" );  // here i use the clone to be fitted because i don t want the main histogram to be modiffied 
	cout << " Chi^2/NDF = " << (double) (FitResult -> Chi2() / FitResult -> Ndf()) << endl;  // here i statistically check how good is the fit spatially.
}
//  =========================================== THE END OF THE ALGORITHM =============================================================================================================================



