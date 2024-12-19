#include "MyLib.hh"
using namespace std;
int a , A;
double miuGen;
std :: random_device RandDev;		// this is the random engine generator declaration
std :: mt19937 Generator(RandDev());
int RandGen( int minn , int maxx )
{
	uniform_int_distribution<> Distr( minn , maxx );
return Distr(Generator);
}
double NeutrMultGenFunc( int a )
{
	if ( (double)a <= 72. )
		return 1.45;
	if ( (double)a > 72. && (double)a <= 74.67 )
		return 36.0576 - 0.489051 * a;
	if ( (double)a > 74.67 && (double)a <= 75.69 )
		return -152.007 + 2.04545 * a;
	if ( (double)a > 75.69 && (double)a <= 77.63 )
		return 76.4675 - 0.974026 * a;
	if ( (double)a > 77.63 && (double)a <= 82.5 )
		return 6.56452 + 0.0782414 * a;
	if ( (double)a > 82.5 && (double)a <= 94. )
		return -6.99 + 0.0857033 * a;
	if ( (double)a > 94. && (double)a <= 98. )
		return 1.04999 + 7.34e-08 * a;
	if ( (double)a > 98. && (double)a <= 102. )
		return -3.89091 + 0.050501 * a;
	if ( (double)a > 102. && (double)a <= 108. )
		return -7.28293 + 0.0834879 * a;
	if ( (double)a > 108. && (double)a <= 123. )
		return -10.3664 + 0.113048 * a;
	if ( (double)a > 123. && (double)a <= 125. )
		return 34.1447 - 0.249827 * a;
	if ( (double)a > 125. && (double)a <= 127. )
		return 108.491 - 0.844156 * a;
	if ( (double)a > 127. && (double)a <= 130. )
		return exp(48.43 - 0.38024 * a);
	if ( (double)a > 130. && (double)a <= 141. )
		return -14.2888 + 0.112928 * a;
	if ( (double)a > 141. && (double)a <= 162. )
		return -5.97558 + 0.0535971 * a;
	if ( (double)a > 162. && (double)a <= 169. )
		return -11.0437 + 0.0850668 * a;
	if ( (double)a > 169. && (double)a <= 174. )
		return -36.7084 + 0.236128 * a;
	if ( (double)a > 174. && (double)a <= 180. )
		return 79.0798 - 0.433105 * a;
	if ( (double)a > 180. )
		return 1.12;
return 0;
}
void NMGen()
{
	int Runs;
	cout << "Give the number of Runs, Runs = ";
	cin >> Runs;
	TH2D *h1 = new TH2D ( "h1" , "Neutron Multiplicity as f(A)" , 500 , 0. , 252. , 50 , 0. , 5. );
	
	while ( Runs )
	{
		A = RandGen( 0 , 252 );
		miuGen = NeutrMultGenFunc( A );
		h1 -> Fill ( (double)A , miuGen );		
		Runs --;
	}
	h1 -> Draw();
}

