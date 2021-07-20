/*

6j.cpp
======

Author:     Grant Bradley
Date:       20 July 2021
Version:    1.0

This file contains the routine 6j adapted from Greg Egan's tenJ routine.


Refs:

*/
#define _USE_MATH_DEFINES
#include <cmath>

//  6j symbol is a real number related to the labelling of edges of a 
//  tetrahedra by SU2 representations. Obtained via combination of four
//  normalized Clebsch-Gordan coefficients along six edges of a tetrahedra
int l_01,l_02,l_03,l_23,l_13,l_12;
int I[4] = {0,1,2,3};
int baseJ1[]={1,1,1,1,1};
I[l_i][l_j] = {l_01,l_02,l_03,
               l_23,l_13,l_12}; //  square this matrix

//  gauge invariant variables given by 6 angles theta_iJ \in [0,\pi]
//  for any triple (I,J,K) of distinct elements == set of all possible
//  spherical tetrahedra 
if (theta_IJ <= theta_IK + theta_JK)
if (2* M_PI >= theta_IJ + theta_IK + theta_JK )




double SIXJfloat sixJ()