#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cfd_types.h"
#include "cfd_discretization.h"
#include "cfd_solvers.h"
#include "cfd_convection_diffusion.h"
#include "cfd_turbulence.h"

static int tr=0,tp=0;
#define T(n) do{tr++;printf("  %s ",n);}while(0)
#define OK do{printf("OK\n");tp++;}while(0)
#define CHK(c) do{if(!(c)){printf("FAIL\n");assert(0);}}while(0)
#define EPS 1e-8

void t1(void){T("grid");CfdGrid1D*g=cfd_grid1d_create(10,1.0);CHK(g);cfd_grid1d_destroy(g);OK;}
void t2(void){T("scalar");CfdScalarField1D*f=cfd_scalar1d_create(5,CFD_FIELD_TEMPERATURE);CHK(f);cfd_scalar1d_fill(f,300);CHK(fabs(f->data[0]-300)<EPS);cfd_scalar1d_destroy(f);OK;}
void t3(void){T("diff");double d=cfd_central_diff_1d(1,3,0.1);CHK(fabs(d-10)<EPS);OK;}
void t4(void){T("2nd");double d=cfd_second_deriv_1d(1,2,3,1);CHK(fabs(d)<EPS);OK;}
void t5(void){T("powerlaw");double aW,aE,aP;cfd_powerlaw_1d(0,1,1,&aW,&aE,&aP);CHK(aW>0);OK;}
void t6(void){T("peclet");double Pe=cfd_cell_peclet(1,0.1,0.05);CHK(fabs(Pe-2)<EPS);OK;}
void t7(void){T("fluid");CfdFluidProperty air;cfd_fluid_air(&air);CHK(air.density>0);OK;}
void t8(void){T("tdma");int n=5;CfdMatrix2D*m=cfd_matrix2d_create(n,1);double dx=0.2;for(int i=0;i<n;i++){m->aW[i]=m->aE[i]=25;m->aP[i]=50;}m->aP[0]=100;m->b[0]=0;m->aP[4]=100;m->b[4]=50;double*x=calloc(n,sizeof(double));CHK(cfd_solve_tdma(m,x,n)==0);free(x);cfd_matrix2d_destroy(m);OK;}
void t9(void){T("convdiff");CfdGrid1D*g=cfd_grid1d_create(20,1);CfdTransportEquation e={1,0.1,0,0,CFD_SCHEME_POWERLAW,0.5};double*p=calloc(20,sizeof(double));CHK(cfd_solve_convdiff_1d(g,&e,0.5,0,1,CFD_SCHEME_POWERLAW,p)==0);free(p);cfd_grid1d_destroy(g);OK;}
void t10(void){T("exact");double p=cfd_convdiff_1d_exact(0.5,1,0,0,1);CHK(fabs(p-0.5)<EPS);OK;}
void t11(void){T("keps");CfdKEpsilonModel ke;cfd_kepsilon_init(&ke);CHK(fabs(ke.C_mu-0.09)<EPS);OK;}
void t12(void){T("wall");double u=cfd_wall_function_uplus(50,0.41,5);CHK(u>10);OK;}
void t13(void){T("dimless");CfdFluidProperty air;cfd_fluid_air(&air);CfdDimensionless d;cfd_fill_dimensionless_full(&air,10,0.1,0.01,0.001,&d);CHK(d.Reynolds>0);OK;}

int main(void){
    printf("CFD Basics Tests\n");
    t1();t2();t3();t4();t5();t6();t7();t8();t9();t10();t11();t12();t13();
    printf("%d/%d passed\n",tp,tr);
    return (tp==tr)?0:1;
}
