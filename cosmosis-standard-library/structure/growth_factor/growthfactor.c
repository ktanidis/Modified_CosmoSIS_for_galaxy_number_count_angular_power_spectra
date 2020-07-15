#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "growthfactor.h"

//Code to calculate the linear growth factor D, and linear growth rate, f. Where a linear perturbation delta grows as
//delta(a') = delta(a)*(D(a')/D(a))^2 and f = dlnD/dlna
//Note anyone using Komatsu's CRL library should note: growth_factor_crl = D *(1+z) and growth_rate_crl = f/(1+z)

double D;
double linear_f;
double w0,wa;
double omega_m,omega_lambda;

#define LIMIT_SIZE 1000

int get_growthfactor(double a,double om, double w, double w2, double *gf)
{
	w0 = w;
	wa = w2;
	omega_m = om;
	omega_lambda = 1.0 - omega_m;
        growth_de(a,gf);
	return 0;
}

// 2 parameter dark energy equation of state

double w (double a)
{
return w0 + (1.0-a)*wa;
}

double w_int(double z, void *param)
{
	return (1. + w(1./(z+1)) )/( 1. + z);
}


double Xde_int (double a,void *params )
{
    if (a == 0.) a = 1e-3;
    double Xde_int = w(a)/a;
    return Xde_int;
}

double Xde (double a)
{
    gsl_integration_workspace * worksp  = gsl_integration_workspace_alloc (LIMIT_SIZE);
    gsl_function F;
    F.function = &Xde_int;
    F.params =0;
    double integral,error;
    if (a == 0.) 
		a = 1e-3;
    gsl_integration_qags (&F, a, 1,  1.0e-20, 1.0e-10, LIMIT_SIZE,worksp, &integral, &error); 
    gsl_integration_workspace_free (worksp);
	
    return omega_m/(1.-omega_m)*exp(-3.*integral);
}

int func (double t, const double y[], double f[], void *params)
{
	double mu = *(double *)params;
	f[0] = y[1];
	f[1] = -( 3.5 - 1.5 * w(t)/( 1 + Xde(t) ) )*y[1]/t - 1.5 * ( 1 - w(t) )/( 1 + Xde(t))*(y[0]/t/t);
	return GSL_SUCCESS;
}
     
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	double mu = *(double *)params;
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2); 
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set (m, 0, 0, 0.0); //dy1/dx1
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 1, 0, - 1.5 * ( 1 - w(t) )/( 1 + Xde(t))*(1./t/t));
	gsl_matrix_set (m, 1, 1, -( 3.5 - 1.5 * w(t)/( 1 + Xde(t) ) )/t);
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}
     
double growth_de (double a, double *gf)
{
	const gsl_odeiv_step_type * T 
		= gsl_odeiv_step_rk4;
     
	gsl_odeiv_step * s 
         = gsl_odeiv_step_alloc (T, 2);
	gsl_odeiv_control * c 
         = gsl_odeiv_control_y_new (1e-6, 0.0);
	gsl_odeiv_evolve * e 
         = gsl_odeiv_evolve_alloc (2);
     
	double mu = 10;
	gsl_odeiv_system sys = {func, jac, 2, &mu};
     
	double t =1.e-3, t1 = a;
	double h = 1e-6;
	double y[2] = {1.,0.}; 
     
	while (t < t1)
	{
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
		if (status != GSL_SUCCESS)
			break;
	}
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

//	return y[0]*a;  //D growth factor
//	return y[1]*a*a/(y[0]*a) +1.; // f = d ln D/ d ln a
	gf[0] = y[0]*a;
	gf[1] = y[1]*a*a/(y[0]*a) +1.;
	return y[0]*a;

}



