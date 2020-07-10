# Modified CosmoSIS for galaxy number count angular power spectra

# Description

The code is a fully-operational modified version of the publicly available package code CosmoSIS ([Zunts et al. 2015](https://arxiv.org/pdf/1409.3409)). In this modular cosmological parameter estimation suite for galaxy number count Limber approximated angular power spectra, we allow for single and multiple tracers, including density fluctuations, redshift-space distortions, and weak lensing magnification.


# Installation

This version of CosmoSIS, similarly to the original one, needs a 64-bit operating system. The modified code is tested only on Ubuntu 16.04 and 18.04. In principle, it can work also on Mac (see <a href="https://bitbucket.org/joezuntz/cosmosis/wiki/Mac%20Installs" target="_blank">**Install on Mac**</a>). After the successful installation the source code should be replaced with this modified version. In addition there are some package dependencies: 

* <a href="https://www.python.org/downloads/release/python-2710/" target="_blank">**python2.7**</a>, or <a href="https://www.python.org/downloads/" target="_blank">**python3.6 or more recent version**</a>

* reasonably recent versions of C/C++/Fortran compilers (GCC 4.8 and above are fine, as are clang/clang++)

* <a href="https://git-scm.com/downloads" target="_blank">**git**</a>

as well as the libraries - <a href="http://mirror.kumi.systems/gnu/gsl/" target="_blank">**gsl 1.16**</a> or above - <a href="https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html" target="_blank">**cfitsio 3.30**</a> or above - <a href="http://www.fftw.org/download.html" target="_blank">**fftw 3**</a> - <a href="http://www.netlib.org/lapack/#_lapack_version_3_5_0" target="_blank">**lapack**</a>. These are accessible with your package manager using apt or yaml

To install python depedencies run: 

```shell 
$ pip install -r config/requirements.txt
```

and in case you encounter a permissions problem add: `$ --user`.

If you wish to run in parallel than you need to install MPI. This is usually pre-installed on clusters and super-computers. Please be sure that you also have mpi4py installed with the same python and MPI:

```shell
$ pip install --no-binary --no-wheel mpi4py
```

Note again that maybe you'll need to add: `$ --user`

After downloding the code, go to the main directory like this:

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

# Framework

The Limber approximated angular power spectrum for galaxy number counts has the following form:

<img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})=\displaystyle \int \frac{W_g ^{A,i}(k_\ell,\chi)W_g ^{B,j}(k_\ell,\chi)}{\chi^2}P_{\text{lin}} (k_\ell)">

with <img src="https://render.githubusercontent.com/render/math?math=\chi"> the comoving diastance at a given redshift <img src="https://render.githubusercontent.com/render/math?math=z">, <img src="https://render.githubusercontent.com/render/math?math=P_{\text{lin}} (k_\ell)"> the linear matter power spectrum at present, <img src="https://render.githubusercontent.com/render/math?math=k_\ell=\frac{\ell %2B 1/2}{\chi}">, <img src="https://render.githubusercontent.com/render/math?math=A,B"> the different tracers and <img src="https://render.githubusercontent.com/render/math?math=i,j"> the indices of the redshift bins. Thus, <img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})"> denotes the correlation between tracer <img src="https://render.githubusercontent.com/render/math?math=A"> in the <img src="https://render.githubusercontent.com/render/math?math=i">th redshift bin and tracer  <img src="https://render.githubusercontent.com/render/math?math=B"> in the <img src="https://render.githubusercontent.com/render/math?math=j">th redshift bin. When  <img src="https://render.githubusercontent.com/render/math?math=A=B">, we have single tracer correlations; similarly, when  <img src="https://render.githubusercontent.com/render/math?math=i=j"> we have same-redshift correlations. The scale and redshift dependent kernel reads:

<img src="https://render.githubusercontent.com/render/math?math=W_{g} ^{A,i}(k_\ell,\chi)=W_{g,den} ^{A,i}(k_\ell,\chi)%2BW_{g,RSD} ^{A,i}(k_\ell,\chi)%2B W_{g,mag} ^{A,i}(k_\ell,\chi)">. The three terms are:

<img src="https://render.githubusercontent.com/render/math?math=W_{g,den} ^{A,i}(k_\ell,\chi)=N_A^i(\chi)b(k_\ell,\chi)D(k_\ell,\chi)">

the galaxy density fluctuations,

<img src="https://render.githubusercontent.com/render/math?math=W_{g,RSD} ^{A,i}(k_\ell,\chi)=\frac{2\ell^2 %2B 2\ell-1}{(2\ell-1)(2\ell%2B 3)}\left[N_A^i(\chi)\right]\Big\{\left[fD\right](k_\ell,\chi)\Big\}-\frac{(\ell-1)\ell}{(2\ell-1)\sqrt{(2\ell-3)(2\ell %2B 1)}}\left[N_A^i (\chi)\right]\Bigg\{\left[ fD\right]\left(k_\ell,\frac{2\ell-3}{2\ell %2B 1}\chi\right)\Bigg\}-\frac{(\ell %2B 1)(\ell %2B 2)}{(2\ell %2B 3)\sqrt{(2\ell %2B 1)(2\ell %2B 5)}}\left[N_A^i(\chi)\right]\Bigg\{\left[fD\right]\left(k_\ell,\frac{2\ell %2B 5}{2\ell %2B 1}\chi\right)\Bigg\}">

the redshift-space distortions, and

<img src="https://render.githubusercontent.com/render/math?math=W_{g,mag} ^{A,i}(k_\ell,\chi)=\frac{3 \Omega_m{H_o}^2}{c^2} \left[1 %2B z(\chi)\right]\chi \widetilde N_A^i(\chi) \left[Q(\chi)-1\right]D(k_\ell,\chi)">

the weak lensing magnification. We denote the linear galaxy bias with <img src="https://render.githubusercontent.com/render/math?math=b(k_\ell,\chi)">; the growth factor with <img src="https://render.githubusercontent.com/render/math?math=D(k_\ell,\chi)"> ; the growth rate of matter perturbations with <img src="https://render.githubusercontent.com/render/math?math=f(k_\ell,\chi)">; and the source redshift distribution of the <img src="https://render.githubusercontent.com/render/math?math=A">-th galaxy population in the <img src="https://render.githubusercontent.com/render/math?math=i">-th redshift bin with <img src="https://render.githubusercontent.com/render/math?math=N_A^i">. Moreover, <img src="https://render.githubusercontent.com/render/math?math=\Omega_m, H_o"> is the total matter fraction in the Universe and the Hubble constant at present respectively, <img src="https://render.githubusercontent.com/render/math?math=c"> the speed of light, while <img src="https://render.githubusercontent.com/render/math?math=Q=5s/2"> is the so-called magnification bias, with <img src="https://render.githubusercontent.com/render/math?math=s"> the slope of the decadic logarithm of the comoving galaxy number density as a function of observed magnitude, taken at the magnitude cut of the survey, and we have defined

<img src="https://render.githubusercontent.com/render/math?math=\tilde N_A^i(\chi)=\displaystyle \int_\chi^\infty \chi^\prime\,\frac{\chi^\prime-\chi}{\chi^\prime}N_A^i(\chi^\prime)">

# Code


# Example

# License
