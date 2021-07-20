/*

theta.cp
========

Author:		Greg Egan
Date:		24 September 2001
Version:	1.0

theta(int twoJ1, int twoJ2, int twoJ3) computes the unnormalised value
of a theta net.

The arguments are three integers, equal to double the spins on the edges of the net.

Reference:	L. Kauffman and S. Lins, Temperley-Lieb Recoupling Theory and
			invariants of 3-Manifolds, Princeton University Press,
			Princeton,  1994.

*/

#include "spin.h"

FACTfloat theta(int twoJ1, int twoJ2, int twoJ3)
{
int sumJ=(twoJ1+twoJ2+twoJ3)/2;
int nn[]={sumJ-twoJ1, sumJ-twoJ2, sumJ-twoJ3, sumJ+1};
int dd[]={twoJ1, twoJ2, twoJ3, 0};
FACTfloat result=multiRatio(nn, dd, 4);
if (sumJ%2==0) return result; else return -result;
}
