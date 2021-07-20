/*

PrimePowers.h
==============

Author:		Greg Egan
Date:		24 September 2001
Version:	1.0

This class represents a rational number as a product of integral powers of primes:

	sign * 2^powers[0] * 3^powers[1] * ... primeList[nPowers-1]^powers[nPowers-1]
	
This enables accurate calculations with ratios of large integers, such as ratios
of factorials.  However, it only permits multiplication and division, since addition
is not efficiently performed on this representation.

The class also manages a list of prime numbers themselves, and a list of
factorisations of factorials.

*/

//	The type "PPfloat" should be defined as the floating point type to which
//	you wish PrimePowers objects to be evaluated.  This would normally be
//	either "double" or "long double".

typedef double PPfloat;
//typedef long double PPfloat;

#pragma longlong on

class PrimePowers
{
public:

//	Data specifying a signed product of integral powers of primes

int sign;					//	-1, 0 or 1
int nPowers;				//	Number of powers in list
int *powers;				//	Array of powers of primes

public:

//	Constructors

PrimePowers();					//	Construct a PrimePowers object with value 1
PrimePowers(long long n);		//	Construct a PrimePowers object with value n
PrimePowers(int sign,			//	Construct a PrimePowers object from caller-supplied
	int nPowers, int *powers);	//	raw data for the sign and powers

//	Destructor

~PrimePowers();

//	Other function declarations

PrimePowers *mult(PrimePowers *inp);	//	Return product of this and inp
PrimePowers *div(PrimePowers *inp);		//	Return quotient of this and inp
PrimePowers *multByFactorials(			//	Return product of this and some
	int *factorials, int *fpowers,		//	integer powers of nfact factorials:
	int nfact);							//	factorials[i]^fpowers[i]
PPfloat evaluate();						//	Value of this, as type PPfloat
long long evaluateLongLong();			//	Value of this, as a long long
PPfloat evaluateSqrt();					//	Value of square root of this, as type PPfloat

long long managedPrimes(int n);			//	The nth prime number, as a long long
PrimePowers *factorialPrimes(int f);	//	PrimePowers object for f!
int *factorise(long long n, int &nf);
};
