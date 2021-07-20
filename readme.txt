readme.txt
==========

Author:		Greg Egan
Date:		24 September 2001
Version:	1.0

This file describes the contents of the "tenJ" package, a set of C++
routines for computing the normalised Riemannian 10j symbols that are
a vital ingredient of the Barrett-Crane model of Riemannian quantum
gravity.  The algorithm used by this package is described in the paper:

	"An Efficient Algorithm for the Riemannian 10j Symbols"
	by J. Daniel Christensen and Greg Egan.

Files in the package
====================

readme.txt					This documentation file.

test.cp						A main program for testing the package.

spin.h						A header file in which some options
							for the package can be set.

factorial.cp				Computes/caches factorials.
theta.cp					Computes theta networks.
tet.cp						Computes tet networks.
tenJ.cp						Computes 10j symbols.

PrimePowers.h				Class definition for the PrimePowers
							class which handles rational numbers
							as a product of integer powers of primes.
PrimePowers.cp				Routines for the PrimePowers class.


Installation
============

Place all the source files in a single directory, compile the six *.cp
files with a C++ compiler of your choice, and link them together.
When you run the program you have created, it should produce output
something like this:

	Mon Sep 24 11:31:31 2001
	Base 2j values:  {{1,1,1,1,1},{1,1,1,1,1}}
	{
	{1, 0.3888888888888888395456433},
	{2, 0.2046666666666666911655881},
	{3, 0.1260272108843536742472935},
	{4, 0.08528970953460747461694069},
	{5, 0.06147982430542724141542266},
	{6, 0.04637444797220449665964281},
	{7, 0.03620568551423208186745839},
	{8, 0.02904438870638861164286126},
	{9, 0.02381753323248704534709219},
	{10, 0.01988902316600982614347437}
	}
	Mon Sep 24 11:31:32 2001
	Elapsed time=1 seconds
	---------------------------------------
	
If it is convenient, you could collect all the object code from the
subroutine files (i.e. all the *.cp files except the test program,
test.cp) into a library.

The routines in each file, and their calling conventions, are
documented in the source code.

Options
=======

These routines should run on any hardware and under any operating
system without alteration, but there are a number of options you can set.
These are mainly to decide what floating point type you wish to use in
various parts of the calculations; higher precision types (such as long
double, if your compiler implements this) can improve the accuracy of
the calculations, especially for very high spins.

(1)	In spin.h

	(a)	USE_PRIME_POWERS		default: true
	
	This can be defined as either true or false, to determine whether
	ratios of factorials are computed with the PrimePowers class, or
	by means of floating point calculations.
	
	(b) FACTfloat				default: double
	
	This is a type definition to decide what floating point type to use
	in factorial routines.  It is also used in the theta routine.
	
	(c) TETfloat				default: double
	
	This is a type definition to decide what floating point type to use
	in the tet routines.
	
	(d) TENJfloat				default: double
	
	This is a type definition to decide what floating point type to use
	in the tenJ routines.
	
	(e) MERGE_TET_THETA			default: true
	
	This can be defined as either true or false, to determine whether
	the ratios of tets and two thetas in the 10j symbol calculations are
	computed separately, or merged into a single calculation.

(2)	In PrimePowers.h

	(a) PPfloat					default: double
	
	This is a type definition to decide what floating point type to use
	when evaluating PrimePowers objects as floating point numbers.
	
(3)	In tenJ.cp

	(a)	logit					default: 0
	
	This can be defined to any positive integer to give a progress message
	after logit steps through the "m1, m2" loops of the 10j calculation.
	This can be useful if you are computing 10j symbols with very large spins,
	which take a while.  A value of 0 gives no progress messages.
	
	(b)	maxC					default: 100
	
	This specifies the maximum dimension of preallocated coefficient matrices.
	For large spins, this parameter might need to be increased; the program will
	print a message to that effect if this is necessary.