/*

PrimePowers.cp
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

#include "PrimePowers.h"
#include <math.h>

#define max(a,b) ((a)>(b))?(a):(b)

//	Static data for a managed list of prime numbers

static long long *primeList=0;			/*	Array of primes	*/
static PPfloat *primeSqrts=0;			/*	Array of square roots of primes */
static int nPrimeList=0;				/*	Count of primes			*/
static int sPrimeList=0;				/*	Size of array primeList	*/
static const int incPrimeList=10000;	/*	Increment when primeList is increased	*/

//	Static data for a managed list of factorised factorials

static PrimePowers **FPlist=0;		/*	Array of PrimePowers describing factorials	*/
static PrimePowers **IPlist=0;		/*	Array of PrimePowers describing integers	*/
static int nFPlist=0;				/*	Count of data in arrays FPlist, IPlist */
static int sFPlist=0;				/*	Size of arrays FPlist, IPlist	*/
static const int incFPlist=100;		/*	Increment when FPlist, IPlist are increased	*/

//	Constructors
//	============

//	Construct a PrimePowers object that represents 1.

PrimePowers::PrimePowers()
{
sign=1;
nPowers=0;
powers=0;
}

//	Construct a PrimePowers object that represents an integer n.

PrimePowers::PrimePowers(long long n)
{
nPowers=0;
powers=0;

if (n==0)
	{
	sign=0;
	return;
	};

if (n<0)				//	set sign
	{
	sign=-1;
	n=-n;
	}
else sign=1;
if (n==1) return;		//	no factors, we're done

if (n<nFPlist)			//	answer is on hand already
	{
	PrimePowers *npp=IPlist[(int)n];
	nPowers=npp->nPowers;
	powers=new int[nPowers];
	for (int i=0;i<nPowers;i++) powers[i]=(npp->powers)[i];
	}
else					//	factorise the integer
	{
	powers=factorise(n, nPowers);
	};
}

//	Construct a PrimePowers object directly from an array

PrimePowers::PrimePowers(int sign, int nPowers, int *powers)
{
this->sign=sign;
this->nPowers=nPowers;
this->powers=powers;
}

//	Destructor
//	==========

PrimePowers::~PrimePowers()
{
if (nPowers!=0 && powers!=0) delete [] powers;
}


//	mult()
//	======

//	Multiply this PrimePowers object by another, and return the result as a PrimePowers
//	object.

PrimePowers *PrimePowers::mult(PrimePowers *inp)
{
int outSign;
if ((outSign=sign*inp->sign)==0) return new PrimePowers(0L);

int inpNPowers=inp->nPowers, outN=max(nPowers, inpNPowers);

//	Create array for result

int *outPowers=new int[outN];

/*	Copy own powers, if any, and zero the rest */

for (int i=0;i<nPowers;i++) outPowers[i]=powers[i];
for (int i=nPowers;i<outN;i++) outPowers[i]=0;

/*	Add powers from other */

int *inpPowers=inp->powers;
for (int i=0;i<inpNPowers;i++) outPowers[i]+=inpPowers[i];

return new PrimePowers(outSign, outN, outPowers);
}

//	div()
//	=====

//	Divide this PrimePowers object by another, and return the result as a PrimePowers
//	object.

PrimePowers *PrimePowers::div(PrimePowers *inp)
{
int outSign;
if ((outSign=sign/inp->sign)==0) return new PrimePowers(0L);

int inpNPowers=inp->nPowers, outN=max(nPowers, inpNPowers);

//	Create array for result

int *outPowers=new int[outN];

/*	Copy own powers, if any, and zero the rest */

for (int i=0;i<nPowers;i++) outPowers[i]=powers[i];
for (int i=nPowers;i<outN;i++) outPowers[i]=0;

/*	Subtract powers from other */

int *inpPowers=inp->powers;
for (int i=0;i<inpNPowers;i++) outPowers[i]-=inpPowers[i];

return new PrimePowers(outSign, outN, outPowers);
}

//	multByFactorials()
//	==================
//
//	Multiply this PrimePowers object by a list of powers of factorials, and return
//	the resulting PrimePowers object.
//
//	factorials[]	lists integers, powers of whose factorials are used to multiply
//
//	fpowers[]		lists the powers for each listed factorial; these can be -ve
//
//	nfact			is the number of factorials listed
//
//	The array returned represents:
//
//	(value of this) * (factorials[0]!)^powers[0] * (factorials[1]!)^powers[1] * ...

PrimePowers *PrimePowers::multByFactorials(int *factorials, int *fpowers, int nfact)
{
if (sign==0) return new PrimePowers(0L);

int outN=nPowers;
for (int i=0;i<nfact;i++) outN=max(factorialPrimes(factorials[i])->nPowers,outN);

//	Create array for result

int *outPowers=new int[outN];

/*	Copy own powers, if any, and zero the rest */

for (int i=0;i<nPowers;i++) outPowers[i]=powers[i];
for (int i=nPowers;i<outN;i++) outPowers[i]=0;

/*	Apply powers from list */

for (int i=0;i<nfact;i++)
	{
	PrimePowers *fpp=factorialPrimes(factorials[i]);
	int n=fpp->nPowers, *a=fpp->powers;
	int fp=fpowers[i];
	for (int j=0;j<n;j++) outPowers[j]+=fp*a[j];
	};

return new PrimePowers(sign, outN, outPowers);
}

/*

evaluate()
==========

Returns the value of this PrimePowers object as a PPfloat.
	
*/

PPfloat PrimePowers::evaluate()
{
int i,j,n;
PPfloat prod=sign, factor;

for (i=0;i<nPowers;i++)
if ((n=powers[i])!=0)
	{
	factor=primeList[i];
	if (n>0) for (j=0;j<n;j++) prod=prod*factor;
	else for (j=0;j<-n;j++) prod=prod/factor;
	};
return prod;
}

/*

evaluateLongLong()
==================

Returns the value of this PrimePowers object as a long long;
this is subject to truncation if there are any negative powers.
	
*/

long long PrimePowers::evaluateLongLong()
{
int i,j,n;
long long prod=sign, factor;

for (i=0;i<nPowers;i++)
if ((n=powers[i])!=0)
	{
	factor=primeList[i];
	if (n>0) for (j=0;j<n;j++) prod=prod*factor;
	else for (j=0;j<-n;j++) prod=prod/factor;
	};
return prod;
}

/*

evaluateSqrt()
==============

Returns the square root of this PrimePowers object, or minus the square root if the
object is negative.

*/

PPfloat PrimePowers::evaluateSqrt()
{
int i,j,n,n2;
PPfloat prod=sign, factor;

for (i=0;i<nPowers;i++)
if ((n=powers[i])!=0)
	{
	factor=primeList[i];
	n2=n/2;
	if (n>0)
		{
		for (j=0;j<n2;j++) prod=prod*factor;
		if (n%2!=0) prod=prod*primeSqrts[i];
		}
	else
		{
		for (j=0;j<-n2;j++) prod=prod/factor;
		if (n%2!=0) prod=prod/primeSqrts[i];
		};
	};
return prod;
}

//	managedPrimes()
//	===============
//
//	Returns the nth entry in a list of prime numbers, creating it if it does
//	not yet exist.

long long PrimePowers::managedPrimes(int n)
{
long long topPrime, c;
int j;
bool isPrime;

if (nPrimeList<=n)
	{
	/*	Work out first candidate to test for being prime	*/
	
	topPrime = nPrimeList==0 ? -1 : primeList[nPrimeList-1];
	if (topPrime>=3) c=topPrime+2;
	else if (topPrime==2) c=3;
	else c=2;
	
	/*	Loop until we get primeList[n]	*/
	
	while (nPrimeList<=n)
		{
		isPrime=true;
		for (j=0;j<nPrimeList;j++)		/*	See if any prime divides c	*/
			{
			long long pj=primeList[j];
			if (pj*pj > c) break;
			if (c % pj == 0)
				{
				isPrime=false;
				break;
				};
			};
			
		/*	We have a new prime to add to the list	*/
			
		if (isPrime)
			{
			
			/*	Filled arrays, so expand them	*/
			
			if (nPrimeList==sPrimeList)
				{
				sPrimeList+=incPrimeList;
				long long *tmp=new long long[sPrimeList];
				PPfloat *td=new PPfloat[sPrimeList];
				if (nPrimeList!=0)
					{
					for (j=0;j<nPrimeList;j++)
						{
						tmp[j]=primeList[j];
						td[j]=primeSqrts[j];
						};
					delete [] primeList;
					delete [] primeSqrts;
					};
				primeList=tmp;
				primeSqrts=td;
				};
			
			primeList[nPrimeList]=c;
			primeSqrts[nPrimeList]=sqrt(c);
			nPrimeList++;
			};
			
		/*	Next candidate to test; normally increment from odd to next odd	*/

		if (c==2) c++; else c+=2;
		};
	};
return primeList[n];
}

//	factorialPrimes()
//	=================
//
//	Returns a pointer to the PrimePowers object for f!, creating it
//	if it does not yet exist.

PrimePowers *PrimePowers::factorialPrimes(int f)
{
int i, j, newSFP;
PrimePowers **tmp;

if (f>=nFPlist)						/*	Don't yet have data for f!	*/
	{
	if (f>=sFPlist)					/*	Make room for lists to go up to f	*/
		{
		i=1+(f+1-sFPlist)/incFPlist;
		newSFP=sFPlist+i*incFPlist;
		
		tmp=new PrimePowers *[newSFP];
		if (FPlist!=0)
			{
			for (j=0;j<sFPlist;j++) tmp[j]=FPlist[j];
			delete [] FPlist;
			};
		FPlist=tmp;
		
		tmp=new PrimePowers *[newSFP];
		if (IPlist!=0)
			{
			for (j=0;j<sFPlist;j++) tmp[j]=IPlist[j];
			delete [] IPlist;
			};
		IPlist=tmp;

		sFPlist=newSFP;
		};
		
	/*	Seed things with 0!, first time we run	*/
	
	if (nFPlist<=0)
		{
		FPlist[0]=new PrimePowers();
		IPlist[0]=new PrimePowers();
		nFPlist=1;
		};
		
	/*	Loop, creating new factorials */
	
	int topFP=nFPlist-1;
	while (topFP<f)
		{
		topFP++;
		int np, *p=factorise(topFP, np);

		IPlist[topFP]=new PrimePowers(1,np,p);
		FPlist[topFP]=IPlist[topFP]->mult(FPlist[topFP-1]);
		};
	nFPlist=topFP+1;
	};
return FPlist[f];
}

//	factorise()
//	===========
//
//	Factorise a long long integer n;
//	returns a pointer to an array which counts powers of prime factors;
//	nf is set equal to the number of entries in the array.

int *PrimePowers::factorise(long long n, int &nf)
{
const int sizeInc=100;
int size=10, *factors=new int[size];
long long d=n;
int mf=-1;
for (int j=0;d>1;j++)
	{
	if (j>=size)
		{
		int newSize=j+sizeInc;
		int *tmp=new int[newSize];
		for (int k=0;k<size;k++) tmp[k]=factors[k];
		size=newSize;
		delete [] factors;
		factors=tmp;
		};
	factors[j]=0;
	long long p=managedPrimes(j);
	while (d % p == 0)
		{
		d = d / p;
		factors[j]++;
		mf=j;
		};
	};
nf=mf+1;
return factors;
}
