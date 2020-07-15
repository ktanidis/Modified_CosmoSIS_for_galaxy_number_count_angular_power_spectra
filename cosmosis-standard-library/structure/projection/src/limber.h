#ifndef _H_LIMBER
#define _H_LIMBER

#include <gsl/gsl_spline.h>
#include <stdbool.h>
#include "interp2d.h"
//#include <gsl/gsl_interp2d.h>
//#include <gsl/gsl_spline2d.h>


// The integrated C_ell are in general allowed to be zero or negative if
// they describe cross-correlations. We use these statuses to describe
// errors where a log was requested also.
// LIMBER_STATUS_NEGATIVE is probably always an error
// but LIMBER_STATUS_ZERO is not necessarily
#define LIMBER_STATUS_OK 0
#define LIMBER_STATUS_ZERO 1
#define LIMBER_STATUS_NEGATIVE 2
#define LIMBER_STATUS_ERROR 3

// These are the options you can set for
// the Limber integrator.
typedef struct limber_config{
	bool xlog;  // The output spline will be in terms of log(ell) not ell
	bool ylog; 
	int n_ell;  // Number of ell values you want in the spline
	double * ell;  // The chosen ell values you want
	double prefactor; //Scaling prefactor
    int status; // did everything go okay?
    double absolute_tolerance;
    double relative_tolerance;
} limber_config;

gsl_spline * limber_integral(limber_config * config, 
	gsl_spline * WbX, gsl_spline * WfX, gsl_spline * WmX, gsl_spline * WbY, gsl_spline * WfY, gsl_spline * WmY, Interpolator2D * P, Interpolator2D * f, Interpolator2D * D, Interpolator2D * BB);

#endif
