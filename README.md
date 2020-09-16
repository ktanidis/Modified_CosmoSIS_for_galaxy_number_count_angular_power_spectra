# Modified CosmoSIS for galaxy number count angular power spectra

# Description

The code is a fully-operational modified version of the publicly available package code CosmoSIS ([Zunts et al. 2015](https://arxiv.org/pdf/1409.3409)). In this modular cosmological parameter estimation suite for galaxy number count Limber approximated angular power spectra, we allow for single and multiple tracers, including density fluctuations, redshift-space distortions, and weak lensing magnification.


# Installation

This version of CosmoSIS, similarly to the original one, needs a 64-bit operating system. The modified code is tested on Ubuntu 16.04 and 18.04 but in principle it should work on all Linux. It can also work on Mac (If you wish to install on Mac, see <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/Mac%20Installs" target="_blank">**Install on Mac**</a>). After the successful installation on Mac the main directory should be replaced with this modified version. 

CosmosSIS needs some package dependencies: 

* <a href="https://www.python.org/downloads/release/python-2710/" target="_blank">**python2.7**</a>, or <a href="https://www.python.org/downloads/" target="_blank">**python3.6 or a more recent version**</a>

* reasonably recent versions of C/C++/Fortran compilers (GCC 4.8 and above are fine, as are clang/clang++)

* <a href="https://git-scm.com/downloads" target="_blank">**git**</a>

as well as the libraries - <a href="http://mirror.kumi.systems/gnu/gsl/" target="_blank">**gsl 1.16**</a> or above - <a href="https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html" target="_blank">**cfitsio 3.30**</a> or above - <a href="http://www.fftw.org/download.html" target="_blank">**fftw 3**</a> and - <a href="http://www.netlib.org/lapack/#_lapack_version_3_5_0" target="_blank">**lapack**</a>. These are accessible with your package manager using apt or yam

To install python depedencies run: 

```shell 
$ pip install -r config/requirements.txt
```

and in case you encounter a permissions problem add: `$ --user`.

If you wish to run in parallel then you need to install MPI. This is usually pre-installed on clusters and super-computers. Please be sure that you also have mpi4py installed with the same python and MPI:

```shell
$ pip install --no-binary --no-wheel mpi4py
```

Note again that maybe you'll need to add: `$ --user`

After downloding the code, go to the main directory like this (optional: you can rename the main directory to 'cosmosis' or whatever name you like):

```shell
$ cd cosmosis
$ gedit setup-my-cosmosis
```

Configure the file `setup-my-cosmosis` by specifying the directories of the required dependencies. Then you should be able to run (inside the `cosmosis` directory):

```shell
$ source setup-my-cosmosis
$ make
```

should you get errors during this step, write us an issue by providing us the `setup-my-cosmosis` file along with the complete output of `make`.

If you quit the terminal you should run again the step:

```shell 
$ source setup-my-cosmosis
```

Check that the installation is successful by running <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/demos/Demo1" target="_blank">**demo1**</a> 

NOTE: In case you encounter a permission denial in a specific directory when you run CosmoSIS, open a new terminal and try `$ chmod +x /directory`. Then you should be able to run the code.

# Framework

The Limber approximated angular power spectrum for galaxy number counts has the following form:

<img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})=\displaystyle \int d\chi \frac{W_g ^{A,i}(k_\ell,\chi)W_g ^{B,j}(k_\ell,\chi)}{\chi^2}P_{\text{lin}} (k_\ell)">

with <img src="https://render.githubusercontent.com/render/math?math=\chi"> the comoving diastance at a given redshift <img src="https://render.githubusercontent.com/render/math?math=z">, <img src="https://render.githubusercontent.com/render/math?math=P_{\text{lin}} (k_\ell)"> the linear matter power spectrum at present, <img src="https://render.githubusercontent.com/render/math?math=k_\ell=\frac{\ell %2B 1/2}{\chi}">, <img src="https://render.githubusercontent.com/render/math?math=A,B"> the different tracers and <img src="https://render.githubusercontent.com/render/math?math=i,j"> the indices of the redshift bins. Thus, <img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})"> denotes the correlation between tracer <img src="https://render.githubusercontent.com/render/math?math=A"> in the <img src="https://render.githubusercontent.com/render/math?math=i">th redshift bin and tracer  <img src="https://render.githubusercontent.com/render/math?math=B"> in the <img src="https://render.githubusercontent.com/render/math?math=j">th redshift bin. When  <img src="https://render.githubusercontent.com/render/math?math=A=B">, we have single tracer correlations; similarly, when  <img src="https://render.githubusercontent.com/render/math?math=i=j"> we have same-redshift correlations. The scale and redshift dependent kernel reads:

<img src="https://render.githubusercontent.com/render/math?math=W_{g} ^{A,i}(k_\ell,\chi)=W_{g,den} ^{A,i}(k_\ell,\chi)%2BW_{g,RSD} ^{A,i}(k_\ell,\chi)%2B W_{g,mag} ^{A,i}(k_\ell,\chi)">. The three terms are:

<img src="https://render.githubusercontent.com/render/math?math=W_{g,den} ^{A,i}(k_\ell,\chi)=N_A^i(\chi)b_A^i(k_\ell,\chi)D(k_\ell,\chi)">

the galaxy density fluctuations,

<img src="https://render.githubusercontent.com/render/math?math=W_{g,RSD} ^{A,i}(k_\ell,\chi)=\frac{2\ell^2 %2B 2\ell-1}{(2\ell-1)(2\ell%2B 3)}N_A^i(\chi)\left[fD\right](k_\ell,\chi)-\frac{(\ell-1)\ell}{(2\ell-1)\sqrt{(2\ell-3)(2\ell %2B 1)}}N_A^i\left(\frac{2\ell-3}{2\ell %2B 1}\chi\right)\left[fD\right]\left(k_\ell,\frac{2\ell-3}{2\ell %2B 1}\chi\right)-\frac{(\ell %2B 1)(\ell %2B 2)}{(2\ell %2B 3)\sqrt{(2\ell %2B 1)(2\ell %2B 5)}}N_A^i\left(\frac{2\ell %2B 5}{2\ell %2B 1}\chi\right)\left[fD \right]\left(k_\ell,\frac{2\ell %2B 5}{2\ell %2B 1}\chi\right)">

the redshift-space distortions, and

<img src="https://render.githubusercontent.com/render/math?math=W_{g,mag} ^{A,i}(k_\ell,\chi)=\frac{3 \Omega_m{H_o}^2}{c^2} \left[1 %2B z(\chi)\right]\chi \widetilde N_A^i(\chi) \left[Q(\chi)-1\right]D(k_\ell,\chi)">

the weak lensing magnification. We denote the linear galaxy bias with <img src="https://render.githubusercontent.com/render/math?math=b(k_\ell,\chi)">; the growth factor with <img src="https://render.githubusercontent.com/render/math?math=D(k_\ell,\chi)"> ; the growth rate of matter perturbations with <img src="https://render.githubusercontent.com/render/math?math=f(k_\ell,\chi)">; and the source redshift distribution of the <img src="https://render.githubusercontent.com/render/math?math=A">-th galaxy population in the <img src="https://render.githubusercontent.com/render/math?math=i">-th redshift bin with <img src="https://render.githubusercontent.com/render/math?math=N_A^i">. Moreover, <img src="https://render.githubusercontent.com/render/math?math=\Omega_m, H_o"> is the total matter fraction in the Universe and the Hubble constant at present respectively, <img src="https://render.githubusercontent.com/render/math?math=c"> the speed of light, while <img src="https://render.githubusercontent.com/render/math?math=Q=5s/2"> is the so-called magnification bias, with <img src="https://render.githubusercontent.com/render/math?math=s"> the slope of the decadic logarithm of the comoving galaxy number density as a function of observed magnitude, taken at the magnitude cut of the survey, and we have defined

<img src="https://render.githubusercontent.com/render/math?math=\tilde N_A^i(\chi)=\displaystyle \int_\chi^\infty \chi^\prime\,\frac{\chi^\prime-\chi}{\chi^\prime}N_A^i(\chi^\prime)">

which is called lensing efficiency.

# Code

To account for single or multi-tracer analysis, we use the module `cosmosis/cosmosis-standard-library/number_density/load_nz/` (see <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/default_modules/load_nz_1">**load_nz**</a>). This module reads from a txt file the <img src="https://render.githubusercontent.com/render/math?math=N_A^i(z)"> distributions for the <img src="https://render.githubusercontent.com/render/math?math=i">th bin of the tracer <img src="https://render.githubusercontent.com/render/math?math=A"> with the format: 1st column redshift <img src="https://render.githubusercontent.com/render/math?math=z">, and the rest the <img src="https://render.githubusercontent.com/render/math?math=N_A^i(z)"> bins. For example for two tracers A,B each having 2 bins in the redshift range <img src="https://render.githubusercontent.com/render/math?math=z"> the columns read : z, tracerA:bin1,tracerA:bin2,tracerB:bin1,tracerB:bin2

The modified part of the code is the module under the directory: `cosmosis/cosmosis-standard-library/structure/projection/src/` (for the original CosmoSIS module version see <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/default_modules/project_2d_1.0">**project_2d**</a>), and more specifically:

* `utils.c`: Loads the function <img src="https://render.githubusercontent.com/render/math?math=P_\text{lin}(k_\ell)">, and calculates <img src="https://render.githubusercontent.com/render/math?math=D(k_\ell,\chi)">, <img src="https://render.githubusercontent.com/render/math?math=f(k_\ell,\chi)">, <img src="https://render.githubusercontent.com/render/math?math=b(k_\ell,\chi)"> (based on eq.5.6 of [Castorina et al. 2014](https://arxiv.org/pdf/1311.1212))

* `kernel.c`: Specifies the considered galaxy number count contributions under the names `DEN` for the galaxy density field, `RSD` for redshift-space distortions and `MAG` for the weak lensing magnification. It also calculates the corresponding normalized <img src="https://render.githubusercontent.com/render/math?math=N_A ^i"> and some prefacors (in the case of weak lensing magnification)

* `limber.c`: Calculates the <img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})"> assuming <img src="https://render.githubusercontent.com/render/math?math=W_{g} ^{A,i}(k_\ell,\chi)=W_{g,den} ^{A,i}(k_\ell,\chi)%2BW_{g,RSD} ^{A,i}(k_\ell,\chi)%2B W_{g,mag} ^{A,i}(k_\ell,\chi)">.

In addition to these, the Python interface of the code is modified as well (the directory `cosmosis/cosmosis-standard-library/structure/projection/`).

* `limber.py`: loads the source code functions

* `project_2d.py`: provides the output <img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})">. Three kernels are used with the names `W_source, F_source` and `M_source` accounting for `DEN, RSD` and `MAG`

NOTE: The current modified CosmoSIS version is valid for galaxy clustering ONLY under the entry `galcl-galcl=source-source-source` (see the Example below). DO NOT attempt to ask output for the section names like CMB_kappa, Shear or Intrinsic alignments, etc.. (for these see again the original module <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/default_modules/project_2d_1.0">**project_2d**</a>).

Finally, we modified the Gaussian likelihood module `cosmosis/cosmosis-standard-library/likelihood/2pt/` to account for the output name `galcl`


# Example

Run the following example `ini` file with : `$ cosmosis examples/my_example/LRGnELG.ini`. This example provides the mock data (data vector and covariance matrix) of multi-tracer galaxy number count angular spectra between the Luminus Red Galaxy and the Emission Line Galaxy samples of the Dark Energy Spectroscopic Instrument ([Aghamousa et al. 2016](https://arxiv.org/pdf/1611.00036)) given a fiducial cosmological model.


```ini
[runtime]
sampler = test

[test]
save_dir=examples/my_example/output
fatal_errors=T

[pipeline]
modules = consistency camb sigma8_rescale load_nz galcl_galcl save_simulation
values = examples/my_example/values.ini
likelihoods =
extra_output = 
quiet=T
timing=F
debug=F

[consistency]
file = cosmosis-standard-library/utility/consistency/consistency_interface.py

[camb]
file = cosmosis-standard-library/boltzmann/camb/camb.so
mode=all
lmax=2500
feedback=0
zmax=1.75
nz=176

[sigma8_rescale]
file = cosmosis-standard-library/utility/sample_sigma8/sigma8_rescale.py

[load_nz]
file = cosmosis-standard-library/number_density/load_nz/load_nz.py
filepath = cosmosis-standard-library/LRGandELG.txt
output_section=nz_source

[galcl_galcl]
file = cosmosis-standard-library/structure/projection/project_2d.py
ell_min = 2.0
ell_max = 1000.0
n_ell = 20
;this is to incorporate density-RSD-magnification for the bins auto and cross-correlations
galcl-galcl = source-source-source
verbose = T

[save_simulation]
file = cosmosis-standard-library/likelihood/2pt/save_2pt.py
galcl_nz_name = source
position_nz_name = source
filename = examples/my_example/LRGnELG.fits
;clobber = T
overwrite = T
make_covariance = T
gaussian_covariance = T
;These values define the survey and the observations being made
;First, some details of the survey itself:
survey_area = 14000.0
number_density_galcl_bin = 0.0288 0.0307 0.0134 0.0022 0.00026243 0.0268 0.0782 0.0725 0.0537 0.1077
number_density_lss_bin = 0.0288 0.0307 0.0134 0.0022 0.00026243 0.0268 0.0782 0.0725 0.0537 0.1077
;we do not care about sigma_e for galaxy clustering
sigma_e_bin = 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
;Then the observations we will generate:
ell_min = 2.0
ell_max = 1000.0
n_ell = 20

data_sets=galcl_cl

```

that reads the input `values.ini` file:

```ini
; Parameters and data in CosmoSIS are organized into sections
; so we can easily see what they mean.
; There is only one section in this case, called cosmological_parameters
[cosmological_parameters]
omega_m      =  0.1 0.3089 0.6
h0           =  0.5 0.6774 1.0
omega_b = 0.0486
sigma8_input = 0.5 0.8159 1.2

;choice for 1 massive (and 2 massless) neutrino with 0.06eV
omnuh2 = 0.00064
massless_nu = 2.046
massive_nu = 1
; Tau (optical depth) is only needed if you need
; the thermal history, or the CMB/matter power
tau = 0.066

; These ones are only needed if you are doing CMB or matter
; power data, not if you just want thermal history or background
; evolution
n_s = 0.9667
A_s = 2.141304e-9
; These parameters can affect everything but are not required - 
; if they are not found in this file they take the sensible
; default values shown here
omega_k = 0.0
w = -1.0
wa = 0.0

;these are the normalisation bias and magnification bias parameters per bin

bias1 = 1.0
bias2 = 1.0
bias3 = 1.0
bias4 = 1.0
bias5 = 1.0
bias6 = 1.0
bias7 = 1.0
bias8 = 1.0
bias9 = 1.0
bias10 = 1.0

alpha_m1 = 1.07
alpha_m2 = 1.10
alpha_m3 = 1.12
alpha_m4 = 1.15
alpha_m5 = 1.17
alpha_m6 = 1.07
alpha_m7 = 1.10
alpha_m8 = 1.12
alpha_m9 = 1.15
alpha_m10 = 1.17

; account for the density fluctuations, RSD and magnification bias contributions. 
;To include the effect specify 1, otherwise 0 (default is 1) 
DEN=1
RSD=1
MAG=1

```
In order to consider or not the density fluctuations, the redshift-space distortions and weak lensing magnification contributions specify the value 1 or 0 respectively for the parameters `DEN, RSD` and `MAG` as seen in the above `values.ini` file. Should you have any inqueries regarding the rest of the CosmoSIS modules and the input cosmological parameters, visit <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/Home">**CosmoSIS_Wiki**</a>.


# License

The code is an available open source and is described in the paper [Tanidis & Camera 2020](https://arxiv.org/pdf/2009.05584). You are welcome to reuse it under the terms of our <a href="https://github.com/ktanidis/Modified_CosmoSIS_for_galaxy_number_count_angular_power_spectra/blob/master/LICENSE">**LICENSE**</a> (<a href="https://opensource.org/licenses/BSD-2-Clause">**BSD-2-Clause**</a>). 

Note that some sub-components of this code like `CAMB`, `CLASS`, `emcee`, `multinest`, etc.. have their own separate licenses that should be taken into account if used. These can be found in the <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/Home">**CosmoSIS_Wiki**</a>.
