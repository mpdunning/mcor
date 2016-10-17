/* file: hppmASub.c
 * Collection of aSub subroutines...
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <time.h>
//#include <sys/time.h>
#include <dbDefs.h>
#include <alarm.h>
#include <registryFunction.h>
#include "aSubRecord.h"
#include <epicsExport.h>
#include <math.h>
#include <stdlib.h>

#define SIZE(x)               (sizeof(x)/sizeof(x[0]))

#define NUMITER 100

typedef unsigned char byte;
typedef unsigned short word;
typedef unsigned int uint;


static double BtoI(double b, double a0, double a1, double a2, double a3, double a4, double a5) {
    double i;
    i=a0 + a1*b +  a2*b*b + a3*b*b*b + a4*b*b*b*b + a5*b*b*b*b*b;
    return i;
}

static long asItoB(aSubRecord *pr) {
/*------------------------------------------------------------------------------
Calculates B(I) by finding a root I(B) using the bisection method.
 *----------------------------------------------------------------------------*/  
    double *i = (double*)pr->a;
    short  *ncoef = (short*)pr->b;
    double *ivbu0 = (double*)pr->c;
    double *ivbu1 = (double*)pr->d;
    double *ivbu2 = (double*)pr->e;
    double *ivbu3 = (double*)pr->f;
    double *ivbu4 = (double*)pr->g;
    double *ivbu5 = (double*)pr->h;
    double *tol = (double*)pr->i;
    double *bmax = (double*)pr->j;
    double *b = (double*)pr->vala;
    
    double b0=0.0,b1=*bmax,ibmax,ib0,ii; 
    int cnt=0;
    
    ibmax=BtoI(*bmax,*ivbu0,*ivbu1,*ivbu2,*ivbu3,*ivbu4,*ivbu5);
    ib0=BtoI(b0,*ivbu0,*ivbu1,*ivbu2,*ivbu3,*ivbu4,*ivbu5);
    
    if(*ncoef<1) {
        printf( "ncoef=%d, must be >= 1\n",*ncoef);
        *b=0.0;
    }
    else if(*ncoef <= 2 && *ivbu1 != 0.0) {
        *b=(*i-*ivbu0)/(*ivbu1);
        //printf("ncoef=%d,i=%f,b=%f\n",*ncoef,*i,*b);
    }
    else {
        while((cnt++) < NUMITER) {
            *b=(b0+b1)/2;
            ii=BtoI(*b,*ivbu0,*ivbu1,*ivbu2,*ivbu3,*ivbu4,*ivbu5) - *i;
            if((ii==0) | ((b1-b0)/2 < *tol)) break;
            //printf("ncoef=%d,cnt=%d,i=%f,b=%f,b0=%f,b1=%f,ii=%f\n",*ncoef,cnt,*i,*b,b0,b1,ii);
            if(ii*ib0 > 0) b0=*b;
            else b1=*b;
        }
    }

    return(0);
}

/**************************************************************************/

epicsRegisterFunction(asItoB);

