/*

tenJ.cp
=======

Author:		Greg Egan
Date:		24 September 2001
Version:	1.0

This file contains two versions of the routine tenJ(); both compute the normalised
value of a 10j symbol, but their different argument signatures distinguish between
versions that accept 10 spins, for the general case, or 1 spin, for the regular case.
In either case, the arguments are integers, which are twice the half-integer spins.

tenJ(int *twoJ1, int *twoJ2) computes the 10j symbol in the general case

tenJ(int twoJ) computes the 10j symbol in the regular case

Reference:	J.D. Christensen and G. Egan, "An Efficient Algorithm for the Riemannian
			10j Symbols".

*/

#include "spin.h"

//	Set logit to the number of inner loop passes at which to print a progress line;
//	set to zero for no progress.

#define logit 0

//	Maximum dimensions for coefficient matrices

#define maxC 100

//	tenJ() computes the normalised value of a 10j symbol
//
//	twoJ1[]		gives double the values of the spins on the five edges joining
//				vertices 0 to 1, 1 to 2, 2 to 3, 3 to 4, 4 to 0
//	twoJ2[]		gives double the values of the spins on the five edges joining
//				vertices 0 to 2, 1 to 3, 2 to 4, 3 to 0, 4 to 1

TENJfloat tenJ(int *twoJ1, int *twoJ2)
{
static int L[5], H[5];					//	Limits on c_i, independent of m's
static int LL[5], HH[5], dim[5];		//	Limits on c_i, taking m's into account
static TENJfloat M[5][maxC][maxC];		//	Coefficient matrices
static TENJfloat v0[maxC], v1[maxC];	//	Vectors used in computing trace

//	Low/high limits on c_i, independent of m's

for (int i=0;i<5;i++)
	{
	int j1=twoJ1[i], j2=twoJ2[i];
	int j1m=twoJ1[mod5(i-1)], j2m2=twoJ2[mod5(i-2)];
	L[i]=max(abs(j1-j2),abs(j1m-j2m2));
	H[i]=min(j1+j2,j1m+j2m2);
	if (1+(H[i]-L[i])/2>maxC)
		{
		printf("Array size maxC=%d exceeded; recompile tenJ.cp with higher maxC\n",maxC);
		exit(1);
		};
	};
	
//	Low/high limits on m's

int mLow, mHigh;
for (int i=0;i<5;i++)
	{
	int j2m=twoJ2[mod5(i-1)], Li=L[i], Hi=H[i];

	int mh=H[i]+j2m;
	if (i==0 || mh<mHigh) mHigh=mh;
	
	int ml;
	if (j2m>=Hi) ml=j2m-Hi;
	else if (j2m<=Li) ml=Li-j2m;
	else ml=(j2m-Li)%2;
	if (i==0 || ml>mLow) mLow=ml;
	};
	
//	Overall sign depends on sum of 2j for all edges, plus the m-independent
//	part of the sign from equation (5) in the paper, (-1)^{2(L_0+j_{2,4})}

int overallParity=L[0]+twoJ2[4];
for (int i=0;i<5;i++)
	{
	overallParity+=(twoJ1[i]+twoJ2[i]);
	};

//	Outermost loop:  for mLow<=m2<=m1<=mHigh

TENJfloat sumOverM=0.0;
int td=(mHigh-mLow)/2+1, tSteps=td*(td+1)/2, steps=0;
for (int m1=mLow;m1<=mHigh;m1+=2)
for (int m2=mLow;m2<=m1;m2+=2)
	{
	if (logit>0 && ++steps%logit==1) printf("Starting step %d of %d\n",steps,tSteps);
	
	//	Low/high limits on c_i, taking current m values into account
	
	bool compatible=true;
	int lowDim=0;		//	Identify index to lowest dimension
	for (int i=0;i<5 && compatible;i++)
		{
		int j2m=twoJ2[mod5(i-1)];
		
		int origLow=L[i], clipLow=max(abs(m1-j2m),abs(m2-j2m));
		int LLi=LL[i]=max(clipLow, origLow);
		
		int origHigh=H[i], clipHigh=min(m1+j2m,m2+j2m);
		int HHi=HH[i]=min(clipHigh, origHigh);
		
		if (HHi<LLi) {compatible=false; break;};
		
		dim[i]=1+(HHi-LLi)/2;
		if (dim[i]<dim[lowDim]) lowDim=i;
		};
	if (!compatible) continue;
	
	//	Compute the M matrices
	
	for (int k=0;k<5;k++)
		{
		int kp1=(k+1)%5;
		int d1=dim[kp1], d2=dim[k];
		int j1=twoJ1[k], j1p=twoJ1[mod5(k+1)];
		int j2=twoJ2[k], j2m=twoJ2[mod5(k-1)], j2p=twoJ2[mod5(k+1)];
		
		for (int i=0;i<d1;i++)
			{
			int ckp=LL[kp1]+2*i;
			#if MERGE_TET_THETA
			TENJfloat factor=(ckp+1);
			for (int j=0;j<d2;j++)
				{
				int ck=LL[k]+2*j;
				M[k][i][j]=
					factor*tetOnThetas(ck,j2,ckp,j2m,m1,j1,j2,ckp,m1,j2m,ckp,j1)*
					tetOnThetas(ck,j2,ckp,j2m,m2,j1,j2,ckp,m2,j2p,ckp,j1p);
				};
			#else
			TENJfloat factor=
				(ckp+1)/
				(theta(j2,ckp,m1)*theta(j2,ckp,m2)*theta(j2m,ckp,j1)*theta(j2p,ckp,j1p));
			for (int j=0;j<d2;j++)
				{
				int ck=LL[k]+2*j;
				M[k][i][j]=
					factor*tet(ck,j2,ckp,j2m,m1,j1)*tet(ck,j2,ckp,j2m,m2,j1);
				};
			#endif
			};
		};
		
	//	Find the trace of their product; outer sum is over basis vectors in the space
	//	with dimensions dim 0 (domain of first matrix, range of last matrix)
	
	TENJfloat trace=0.0;
	int d0=dim[lowDim], d1=dim[mod5(lowDim+1)];
	for (int l0=0;l0<d0;l0++)
		{
		//	Set v0 to product of first matrix and basis vector
		
		for (int i=0;i<d1;i++) v0[i]=M[lowDim][i][l0];
		
		//	Multiply by the next three matrices
		
		TENJfloat *vIn=v0, *vOut=v1, *tmp;
		for (int k=1;k<=3;k++)
			{
			//	Set vout = M[k] vin
			
			int ks=mod5(lowDim+k), dIn=dim[ks];
			int ksp=mod5(lowDim+k+1), dOut=dim[ksp];
			for (int i=0;i<dOut;i++)
				{
				TENJfloat vs=0.0;
				for (int j=0;j<dIn;j++) vs+=M[ks][i][j]*vIn[j];
				vOut[i]=vs;
				};
			
			//	Swap input/output vectors
			
			tmp=vIn;
			vIn=vOut;
			vOut=tmp;
			};
		
		//	Add into trace only the relevant coordinate of product with final matrix
		
		int m4=mod5(lowDim+4), dIn=dim[m4];
		for (int j=0;j<dIn;j++) trace+=M[m4][l0][j]*vIn[j];
		};
	
	//	Accumulate into the sum over the m's
	
	TENJfloat term=(m1+1)*(m2+1)*trace*
		((overallParity-(m1+m2)/2)%2==0?1:-1);
	if (m1!=m2) term*=2;
	sumOverM+=term;
	};

return sumOverM;
}


//	Version for regular 10j symbol

TENJfloat tenJ(int twoJ)
{
int H;									//	Limits on c_i, independent of m's
int LL, HH, dim;						//	Limits on c_i, taking m's into account
static TENJfloat M[maxC][maxC];			//	Coefficient matrices
static TENJfloat v0[maxC], v1[maxC];		//	Vectors used in computing trace

//	Low/high limits on c_i, independent of m's

H=2*twoJ;
if (1+H/2>maxC)
	{
	printf("Array size maxC=%d exceeded; recompile tenJ.cp with higher maxC\n",maxC);
	exit(1);
	};
	
//	Low/high limits on m's

int mLow=twoJ%2, mHigh=3*twoJ;

//	Overall sign depends on sum of 2j for all edges, plus the m-independent
//	part of the sign from equation (5) in the paper, (-1)^{2(L_0+j_{2,4})}.
//
//	In the regular case, L_0 is zero, and all 10 edges have the same 2j value,
//	so the overall sign comes from the parity of 11*twoJ.

int overallParity=11*twoJ;

//	Outermost loop:  for mLow<=m2<=m1<=mHigh

TENJfloat sumOverM=0.0;
int td=(mHigh-mLow)/2+1, tSteps=td*(td+1)/2, steps=0;
for (int m1=mLow;m1<=mHigh;m1+=2)
for (int m2=mLow;m2<=m1;m2+=2)
	{
	if (logit>0 && ++steps%logit==1) printf("Starting step %d of %d\n",steps,tSteps);
	
	//	Low/high limits on c_i, taking current m values into account
	
	LL=max(abs(m1-twoJ),abs(m2-twoJ));
	HH=min(m2+twoJ, H);
	if (HH<LL) continue;
	
	dim=1+(HH-LL)/2;
	
	//	Compute the M matrix
	
	for (int i=0;i<dim;i++)
		{
		int ckp=LL+2*i;
		#if MERGE_TET_THETA
		TENJfloat factor=(ckp+1);
		for (int j=0;j<dim;j++)
			{
			int ck=LL+2*j;
			M[i][j]=
				factor*tetOnThetas(ck,twoJ,ckp,twoJ,m1,twoJ,twoJ,ckp,m1,twoJ,ckp,twoJ)*
				tetOnThetas(ck,twoJ,ckp,twoJ,m2,twoJ,twoJ,ckp,m2,twoJ,ckp,twoJ);
			};
		#else
		TENJfloat factor=
			(ckp+1)/
			(theta(twoJ,ckp,m1)*theta(twoJ,ckp,m2)*theta(twoJ,ckp,twoJ)*theta(twoJ,ckp,twoJ));
		for (int j=0;j<dim;j++)
			{
			int ck=LL+2*j;
			M[i][j]=
				factor*tet(ck,twoJ,ckp,twoJ,m1,twoJ)*tet(ck,twoJ,ckp,twoJ,m2,twoJ);
			};
		#endif
		};
		
	//	Find the trace of their product; outer sum is over basis vectors in the space
	//	with dimensions dim 0 (domain of first matrix, range of last matrix)
	
	TENJfloat trace=0.0;
	for (int l0=0;l0<dim;l0++)
		{
		//	Set v0 to product of first matrix and basis vector
		
		for (int i=0;i<dim;i++) v0[i]=M[i][l0];
		
		//	Multiply by the next three matrices
		
		TENJfloat *vIn=v0, *vOut=v1, *tmp;
		for (int k=1;k<=3;k++)
			{
			//	Set vout = M[k] vin
			
			for (int i=0;i<dim;i++)
				{
				TENJfloat vs=0.0;
				for (int j=0;j<dim;j++) vs+=M[i][j]*vIn[j];
				vOut[i]=vs;
				};
			
			//	Swap input/output vectors
			
			tmp=vIn;
			vIn=vOut;
			vOut=tmp;
			};
		
		//	Add into trace only the relevant coordinate of product with final matrix
		
		for (int j=0;j<dim;j++) trace+=M[l0][j]*vIn[j];
		};
	
	//	Accumulate into the sum over the m's
	
	TENJfloat term=(m1+1)*(m2+1)*trace*
		((overallParity-(m1+m2)/2)%2==0?1:-1);
	if (m1!=m2) term*=2;
	sumOverM+=term;
	};

return sumOverM;
}
