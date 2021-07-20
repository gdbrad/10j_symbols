/*

tet.cp
======

Author:		Greg Egan
Date:		24 September 2001
Version:	1.0

This file contains two routines, tet() and tetOnThetas().  Both return results of
type TETfloat, which is just a standard floating point type, chosen in the header
file "spin.h".

tet(int a, int b, int c, int d, int e, int f) computes the unnormalised value
of a tetrahedral net:

       ______________
       |             |
      / \            |
     /   \           |
  b /     \ c        |
   /   f   \         |
   ---------         | e
   \       /         |
    \     /          |
   a \   / d         |
      \ /            |
       |             |
       ______________


This is written as Tet[a b e]
                      [c d f]

These arguments are all integers, equal to double the values of the spins on
the edges of the net.

tetOnThetas(int a, int b, int c, int d, int e, int f,
	int twoJ1a, int twoJ2a, int twoJ3a, int twoJ1b, int twoJ2b, int twoJ3b)
computes the value of the same unnormalised tetrahedral net, divided by two
theta nets, in a form that is useful for the 10j symbol calculations, and which
offers greater possibilities for the cancellation of factorials.

Reference:	L. Kauffman and S. Lins, Temperley-Lieb Recoupling Theory and
			invariants of 3-Manifolds, Princeton University Press,
			Princeton,  1994.

*/

#include "spin.h"

//	Declarations for private helper functions in this file

int arrayMin(int *array, int size);
int arrayMax(int *array, int size);

//	tet()
//	=====

TETfloat tet(int a, int b, int c, int d, int e, int f)
{
int aa[]={(a+b+f)/2, (b+c+e)/2, (c+d+f)/2, (a+d+e)/2};
int bb[]={(b+d+e+f)/2, (a+c+e+f)/2, (a+b+c+d)/2};

//	Compute ratios between consecutive terms in the sum, and identify when this ratio
//	is closest to -1, which should be at the peak in absolute value

int sumLo=arrayMax(aa,4), sumHi=arrayMin(bb,3), nterms=sumHi-sumLo+1;
TETfloat *ratios=new TETfloat[nterms-1];
int r=0, ls=sumLo;
TETfloat lg=1e30;
for (int s=sumLo+1;s<=sumHi;s++)
	{
	int sm=s-1;
	TETfloat rr=ratios[r++]=
				-(TETfloat)((s+1)*(bb[0]-sm)*(bb[1]-sm)*(bb[2]-sm))
				/(TETfloat)((s-aa[0])*(s-aa[1])*(s-aa[2])*(s-aa[3]));

	TETfloat g=abs(rr+1);
	if (g<lg) {ls=s; lg=g;};
	};

//	Pull out the peak term as part of the common factor, and compute the overall
//	common factor.

//	Factorials in the numerator of the common factor
	
int nn[]={
	bb[0]-aa[0], bb[0]-aa[1], bb[0]-aa[2], bb[0]-aa[3],
	bb[1]-aa[0], bb[1]-aa[1], bb[1]-aa[2], bb[1]-aa[3], 
	bb[2]-aa[0], bb[2]-aa[1], bb[2]-aa[2], bb[2]-aa[3], ls+1};

//	Factorials in the denominator of the common factor

int dd[]={a,b,c,d,e,f,
	ls-aa[0], ls-aa[1], ls-aa[2], ls-aa[3],
	bb[0]-ls, bb[1]-ls, bb[2]-ls
	};

TETfloat commonFactor=multiRatio(nn, dd, 13);
	
//	Set the sign of the common factor

if (ls%2==1) commonFactor=-commonFactor;

//	Sum all terms, with the common factor removed, starting from
//	the term that gave the common factor itself.

TETfloat sum=1.0, term1=1.0, term2=1.0;
int s1=ls+1, s2=ls-1;
while (true)
	{
	bool ok1=(s1<=sumHi), ok2=(s2>=sumLo);
	if (!(ok1||ok2)) break;
	if (ok1)
		{
		term1*=ratios[s1-sumLo-1];
		sum+=term1;
		s1++;
		};
	if (ok2)
		{
		term2/=ratios[s2-sumLo];
		sum+=term2;
		s2--;
		};
	};

delete [] ratios;
TETfloat result=sum*commonFactor;

return result;
}

//	tetOnThetas()
//	=============

TETfloat tetOnThetas(int a, int b, int c, int d, int e, int f,
	int twoJ1a, int twoJ2a, int twoJ3a, int twoJ1b, int twoJ2b, int twoJ3b)
{
int aa[]={(a+b+f)/2, (b+c+e)/2, (c+d+f)/2, (a+d+e)/2};
int bb[]={(b+d+e+f)/2, (a+c+e+f)/2, (a+b+c+d)/2};

//	Compute ratios between consecutive terms in the sum, and identify when this ratio
//	is closest to -1, which should be at the peak in absolute value

int sumLo=arrayMax(aa,4), sumHi=arrayMin(bb,3), nterms=sumHi-sumLo+1;
TETfloat *ratios=new TETfloat[nterms-1];
int r=0, ls=sumLo;
TETfloat lg=1e30;
for (int s=sumLo+1;s<=sumHi;s++)
	{
	int sm=s-1;
	TETfloat rr=ratios[r++]=
				-(TETfloat)((s+1)*(bb[0]-sm)*(bb[1]-sm)*(bb[2]-sm))
				/(TETfloat)((s-aa[0])*(s-aa[1])*(s-aa[2])*(s-aa[3]));

	TETfloat g=abs(rr+1);
	if (g<lg) {ls=s; lg=g;};
	};

//	Pull out the peak term as part of the common factor, and compute the overall
//	common factor, divided by the specified thetas.

int sumJa=(twoJ1a+twoJ2a+twoJ3a)/2;
int sumJb=(twoJ1b+twoJ2b+twoJ3b)/2;

//	Factorials in the numerator of the common factor
	
int nn[]={
	bb[0]-aa[0], bb[0]-aa[1], bb[0]-aa[2], bb[0]-aa[3],
	bb[1]-aa[0], bb[1]-aa[1], bb[1]-aa[2], bb[1]-aa[3], 
	bb[2]-aa[0], bb[2]-aa[1], bb[2]-aa[2], bb[2]-aa[3],
	ls+1,
	twoJ1a, twoJ2a, twoJ3a, 0,
	twoJ1b, twoJ2b, twoJ3b, 0
	};

//	Factorials in the denominator of the common factor

int dd[]={a,b,c,d,e,f,
	ls-aa[0], ls-aa[1], ls-aa[2], ls-aa[3],
	bb[0]-ls, bb[1]-ls, bb[2]-ls,
	sumJa-twoJ1a, sumJa-twoJ2a, sumJa-twoJ3a, sumJa+1,
	sumJb-twoJ1b, sumJb-twoJ2b, sumJb-twoJ3b, sumJb+1
	};

TETfloat commonFactor=multiRatio(nn, dd, 21);
	
//	Set the sign of the common factor

if ((ls+sumJa+sumJb)%2==1) commonFactor=-commonFactor;

//	Sum all terms, with the common factor removed, starting from
//	the term that gave the common factor itself.

TETfloat sum=1.0, term1=1.0, term2=1.0;
int s1=ls+1, s2=ls-1;
while (true)
	{
	bool ok1=(s1<=sumHi), ok2=(s2>=sumLo);
	if (!(ok1||ok2)) break;
	if (ok1)
		{
		term1*=ratios[s1-sumLo-1];
		sum+=term1;
		s1++;
		};
	if (ok2)
		{
		term2/=ratios[s2-sumLo];
		sum+=term2;
		s2--;
		};
	};

delete [] ratios;
TETfloat result=sum*commonFactor;

return result;
}

//	Minimum and maximum values in an integer array

int arrayMin(int *array, int size)
{
int m=array[0];
for (int i=1;i<size;i++)
	{
	int a=array[i];
	if (a<m) m=a;
	};
return m;
}

int arrayMax(int *array, int size)
{
int m=array[0];
for (int i=1;i<size;i++)
	{
	int a=array[i];
	if (a>m) m=a;
	};
return m;
}

