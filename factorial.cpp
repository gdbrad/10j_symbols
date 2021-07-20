/*

factorial.cp
============

Author:		Greg Egan
Date:		24 September 2001
Version:	1.0

This file contains two implementations of a routine:

	FACTfloat multiRatio(int *nn, int *dd, int size)
	
The argument arrays nn[] and dd[] list numbers whose factorials are to be
computed in the numerator and denominator of the final result; "size" gives
the number of entries.  If one array has fewer entries than the other, it should
be padded out with zeroes.  Thus the result returned is:

	nn[0]! nn[1]! ... nn[size-1]!
	-----------------------------
	dd[0]! dd[1]! ... dd[size-1]!
	
The routine attempts to exploit cancellation between the numerator and denominator
to avoid underflow or overflow and maximise accuracy.

The first implementation makes use of "PrimePowers" objects, which represent rational
numbers as integer powers of primes.

The second implementation performs all computations using ordinary floating point
arithmetic, of type FACTfloat.

In both implementations, a cache of factorials is accumulated (in the PrimePowers
implementation, this cache is handled in the PrimePowers code).

*/

#include "spin.h"

//	-----------------------------------
//	**** PrimePowers implementation ***
//	-----------------------------------

#if USE_PRIME_POWERS

#include "PrimePowers.h"

//	multiRatio()
//	============
//
//	Return the product of several ratios of factorials

FACTfloat multiRatio(int *nn, int *dd, int size)
{
//	Convert the lists of factorials in the numerator and denominator
//	into a single list with positive and negative powers.

int *f=new int[2*size];
int *p=new int[2*size];
int n=0;
for (int i=0;i<size;i++)
	{
	if (nn[i]!=0)
		{
		f[n]=nn[i];
		p[n++]=1;
		};
	if (dd[i]!=0)
		{
		f[n]=dd[i];
		p[n++]=-1;
		};
	};
	
//	Create and evaluate a PrimePowers object of the required number.

PrimePowers *a=new PrimePowers();
PrimePowers *b=a->multByFactorials(f,p,n);
FACTfloat result=(FACTfloat)b->evaluate();

delete [] f;
delete [] p;
delete a;
delete b;
return result;
}

#else

//	--------------------------------------
//	**** Floating Point implementation ***
//	--------------------------------------

FACTfloat factorialRatio(int f, int g);

//	factorial() returns factorial of non-negative integers.

FACTfloat factorial(int f)
{
return factorialRatio(f,0);	//	Obtain factorial from cache as f!/0!
}

//	factorialRatio() computes/caches ratios of factorials of non-negative integers.

FACTfloat **fList=0;		/*	Array of factorial ratios	*/
int topF=-1;				/*	Highest factorial for which we have data */
int sFlist=0;				/*	Current size of array fList	*/
int incFlist=100;			/*	Increment in size for fList */

FACTfloat factorialRatio(int f, int g)
{
if (f==g) return 1.0;		/* f!/f! */

//	Ensure f>g, and set flag if we need to invert the result

bool invert=false;
if (g>f)
	{
	int t=f;
	f=g;
	g=t;
	invert=true;
	};

int i, j, nsFlist;
FACTfloat **tmp;

if (f>topF)							/*	Don't yet have data for f!	*/
	{
	if (f>=sFlist)					/*	Make room for lists to go up to f	*/
		{
		i=1+(f+1-sFlist)/incFlist;
		nsFlist=sFlist+i*incFlist;
		
		tmp=new FACTfloat *[nsFlist];
		if (fList!=0)
			{
			for (j=0;j<sFlist;j++) tmp[j]=fList[j];
			delete [] fList;
			};
		fList=tmp;
		
		sFlist=nsFlist;
		};
		
	/*	Seed things with 0!, first time we run	*/
	
	if (topF<0)
		{
		topF=0;
		fList[0]=new FACTfloat[1];
		fList[0][0]=1.0;
		};
		
	/*	Loop, creating new sets of factorial ratios */
	
	while (topF<f)
		{
		topF++;
		FACTfloat *latest=fList[topF]=new FACTfloat[topF+1], *prev=fList[topF-1];
		for (int i=0;i<topF;i++) latest[i]=prev[i]*topF;
		latest[topF]=1.0;
		};
	};
return invert ? 1.0/fList[f][g] : fList[f][g];
}

//	Compare two integers for qsort()

int cmp(const void *a, const void *b)
{
return *((int *)a) - *((int *)b);
}

//	multiRatio()
//	============
//
//	Return the product of several ratios of factorials

FACTfloat multiRatio(int *nn, int *dd, int size)
{
//	Sort numerator and denominator data into order, to improve cancellation.

qsort(nn, size, sizeof(int), cmp);
qsort(dd, size, sizeof(int), cmp);

//	Form the individual ratios of factorials

FACTfloat *ratios=new FACTfloat[size];
bool *used=new bool[size];

int extraN=0, extraD=0, nOK=0;
for (int k=0;k<size;k++)
	{
	ratios[k]=factorialRatio(nn[k],dd[k]);
	
	//	Record info about bad ratios to remedy later
	
	if (ratios[k]==0.0 || isinf(ratios[k]) || isnan(ratios[k]))
		{
		ratios[k]=0.0;			//	Mark ratio as bad
		used[k]=true;			//	Mark as already used, so we don't use
		if (nn[k]>dd[k])		//	Count up extra factors needed in numerator/denom.
			extraN+=nn[k]-dd[k];
		else
			extraD+=dd[k]-nn[k];
		}
	else
		{
		used[k]=false;			//	Mark as not yet used
		nOK++;					//	Count up valid ratios
		};
	};
	
//	Multiply all valid ratios together, avoiding overflow/underflow

FACTfloat result=1.0, p;
for (int k=0;k<nOK;k++)
	{
	//	Scan for an unused ratio to multiply by
	
	for (int l=0;l<size;l++)
		if (!used[l])
			{
			p=result*ratios[l];
			
			//	If product is bad, don't store it
			
			if (p==0.0 || isinf(p) || isnan(p)) continue;
			
			//	Product is OK, so store it and mark ratio as used
			
			result=p;
			used[l]=true;
			break;
			};
	};
	
//	Add any unused ratios to the list of problem ratios
	
for (int k=0;k<size;k++)
	if (!used[k])
		{
		ratios[k]=0.0;				//	Mark ratio as bad
		if (nn[k]>dd[k])			//	Count up extra factors needed in numerator/denom.
			extraN+=nn[k]-dd[k];
		else
			extraD+=dd[k]-nn[k];
		nOK--;
		};

//	Deal with any factorial ratios that couldn't be computed as a whole due to
//	overflow or underflow.

if (nOK<size)
	{
	//	Create arrays to hold any extra factors needed in the numerator and denominator
	
	int *en=extraN==0 ? NULL : new int[extraN];
	int *ed=extraD==0 ? NULL : new int[extraD];
	int nen=0, ned=0;
	
	//	Gather all the extra factors in the numerator and denominator

	for (int k=0;k<size;k++)
	if (ratios[k]==0)
		{
		if (nn[k]>dd[k])
			{
			for (int l=nn[k];l>dd[k];l--) en[nen++]=l;
			}
		else
			{
			for (int l=dd[k];l>nn[k];l--) ed[ned++]=l;
			};
		};
		
	//	Multiply them into the common factor, avoiding underflow or overflow
		
	while (nen>0 || ned>0)
		{
		if (result>1.0 && ned>0) result/=ed[--ned];
		else if (result<1.0 && nen>0) result*=en[--nen];
		else if (ned>0) result/=ed[--ned];
		else if (nen>0) result*=en[--nen];
		};
	
	if (en!=NULL) delete [] en;
	if (ed!=NULL) delete [] ed;
	};

delete [] ratios;
delete [] used;
return result;
}

#endif
