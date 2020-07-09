# Modified CosmoSIS for galaxy number count angular power spectra

# Description

The code is a fully-operational modified version of the publicly available package code CosmoSIS ([Zunts et al. 2015](https://arxiv.org/pdf/1409.3409)). In this modular cosmological parameter estimation suite for galaxy number count Limber approximated angular power spectra, we allow for single and multiple tracers, including density fluctuations, redshift-space distortions, and weak lensing magnification.

# Framework

The Limber approximated angular power spectrum for galaxy number counts has the following form:

<img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})=\bigint \frac{W_g ^{A,i}(k_\ell,\chi)W_g ^{B,j}(k_\ell,\chi)}{\chi^2}P_{lin} (k_\ell)">

with <img src="https://render.githubusercontent.com/render/math?math=P_{lin} (k_\ell)"> the linear matter power spectrum at present, <img src="https://render.githubusercontent.com/render/math?math=k_\ell=\frac{\ell + 1/2}{\chi}">, <img src="https://render.githubusercontent.com/render/math?math=A,B"> the different tracers and <img src="https://render.githubusercontent.com/render/math?math=i,j"> the indices of the redshift bins. Thus, <img src="https://render.githubusercontent.com/render/math?math=C_{\ell\gg 1} ^g ({z_i ^A,z_j ^B})"> denotes the correlation between tracer <img src="https://render.githubusercontent.com/render/math?math=A"> in the <img src="https://render.githubusercontent.com/render/math?math=i">th redshift bin and tracer  <img src="https://render.githubusercontent.com/render/math?math=B"> in the <img src="https://render.githubusercontent.com/render/math?math=j">th redshift bin. When  <img src="https://render.githubusercontent.com/render/math?math=A=B">, we have single tracer correlations; similarly, when  <img src="https://render.githubusercontent.com/render/math?math=i=j"> we have same-redshift correlations.
