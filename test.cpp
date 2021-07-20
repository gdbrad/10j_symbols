/*

test.cp
=======

Author:		Greg Egan
Date:		24 September 2001
Version:	1.0

This file contains a main program that uses the "tenJ" package routines to print
some test values for 10j symbols.

*/

#include "spin.h"

//	File in which to log results; comment out to send results to console only.

//	#define logfile "10j.dat"

//	Choose a set of base spins by commenting out all but one set.

//	baseJ1[]	gives double the values of the spins on the five edges joining
//				vertices 0 to 1, 1 to 2, 2 to 3, 3 to 4, 4 to 0
//	baseJ2[]	gives double the values of the spins on the five edges joining
//				vertices 0 to 2, 1 to 3, 2 to 4, 3 to 0, 4 to 1

int baseJ1[]={1,1,1,1,1}, baseJ2[]={1,1,1,1,1};
//int baseJ1[]={2,2,2,2,2}, baseJ2[]={1,1,1,1,1};
//int baseJ1[]={2,1,1,2,1}, baseJ2[]={1,1,1,2,2};
//int baseJ1[]={2,3,1,2,1}, baseJ2[]={1,3,1,2,2};
//int baseJ1[]={3,1,1,3,1}, baseJ2[]={1,1,1,1,1};

//	Minimum, maximum multiples of base spins, and increment for loop

#define MINMULT 1
#define MAXMULT 10
#define INCMULT 1

//	Output a string both to console and log file;
//	log file is opened and closed fresh each time to ensure no
//	data is lost on a crash.

void output(char *string)
{
printf("%s",string);

#ifdef logfile
	FILE *fp=fopen(logfile,"aa");
	fprintf(fp,"%s",string);
	fclose(fp);
#endif
}

static char buffer[512];

//	Main program
//	============

int main()
{
time_t t1, t2;
struct tm *lt;
int twoJ1[5], twoJ2[5];

t1=time(NULL);
lt=localtime(&t1);
output(asctime(lt));

sprintf(buffer,"Base 2j values:  {{%d,%d,%d,%d,%d},{%d,%d,%d,%d,%d}}\n",
	baseJ1[0],baseJ1[1],baseJ1[2],baseJ1[3],baseJ1[4],
	baseJ2[0],baseJ2[1],baseJ2[2],baseJ2[3],baseJ2[4]);
output(buffer);

bool regular=true;
for (int i=0;i<5;i++) regular&=(baseJ1[i]==baseJ1[0] && baseJ2[i]==baseJ1[0]);

output("{\n");
for (int i=MINMULT;i<=MAXMULT;i+=INCMULT)
	{
	double res;
	if (regular) res=tenJ(i*baseJ1[0]); 
	else
		{
		for (int k=0;k<5;k++)
			{
			twoJ1[k]=i*baseJ1[k];
			twoJ2[k]=i*baseJ2[k];
			};
		res=tenJ(twoJ1,twoJ2);
		};
	sprintf(buffer,"{%d, %.25lg}",i,res);
	output(buffer);
	if (i==MAXMULT) output("\n");
	else output(",\n");
	};
output("}\n");
t2=time(NULL);
lt=localtime(&t2);
output(asctime(lt));
sprintf(buffer,"Elapsed time=%d seconds\n",t2-t1);
output(buffer);
output("---------------------------------------\n");
}
