#ifndef _H_LIMBER_UTILS
#define _H_LIMBER_UTILS

#include "gsl/gsl_spline.h"
#include "cosmosis/datablock/c_datablock.h"
#include "interp2d.h"


gsl_spline * spline_from_arrays(int n, double * x, double *y);

DATABLOCK_STATUS save_spline(c_datablock * block, const char * section, 
	const char * n_name, const char * x_name, const char * y_name,
	gsl_spline * s);

gsl_spline * load_spline(c_datablock * block, const char * section, 
	const char * x_name, const char * y_name);



void reverse(double * x, int n);


Interpolator2D * 
load_interpolator(c_datablock * block,
	const char * section,
	const char * k_name, const char * z_name, const char * P_name);


Interpolator2D * 
load_interpolator_chi(c_datablock * block, gsl_spline * chi_of_z,
	const char * section,
	const char * k_name, const char * z_name, const char * P_name);

Interpolator2D * 
load_interpolator_function(c_datablock * block, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);

Interpolator2D * 
load_interpolator_chi_function(c_datablock * block, 
	gsl_spline * chi_of_z_spline,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);



//-------------------for scale dependent growth rate

Interpolator2D * 
load_interpolator2(c_datablock * block,
	const char * section, 
	const char * k_name, const char * z_name, const char * T_name);


Interpolator2D * 
load_interpolator_chi2(c_datablock * block, gsl_spline * chi_of_z, gsl_spline * zi_of_chi,
	const char * section, 
	const char * k_name, const char * z_name, const char * T_name);

Interpolator2D * 
load_interpolator_function2(c_datablock * block, 
	const char * section, 
	const char * k_name, const char * z_name, const char * T_name,
	interp2d_modifier_function function, void * args);

Interpolator2D * 
load_interpolator_chi_function2(c_datablock * block, 
	gsl_spline * chi_of_z_spline, gsl_spline * zi_of_chi_spline,
	const char * section, 
	const char * k_name, const char * z_name, const char * T_name,
	interp2d_modifier_function function, void * args);

//-------------------for scale dependent growth factor

Interpolator2D * 
load_interpolator3(c_datablock * block,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name);


Interpolator2D * 
load_interpolator_chi3(c_datablock * block, gsl_spline * chi_of_z, gsl_spline * zi_of_chi,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name);

Interpolator2D * 
load_interpolator_function3(c_datablock * block, 
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);

Interpolator2D * 
load_interpolator_chi_function3(c_datablock * block, 
	gsl_spline * chi_of_z_spline, gsl_spline * zi_of_chi_spline,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);



//-------------------for scale dependent galaxy bias

Interpolator2D * 
load_interpolator4(c_datablock * block,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name);


Interpolator2D * 
load_interpolator_chi4(c_datablock * block, gsl_spline * chi_of_z, gsl_spline * zi_of_chi,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name);

Interpolator2D * 
load_interpolator_function4(c_datablock * block, 
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);

Interpolator2D * 
load_interpolator_chi_function4(c_datablock * block, 
	gsl_spline * chi_of_z_spline, gsl_spline * zi_of_chi_spline,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);


#endif
