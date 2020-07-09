# Modified CosmoSIS for galaxy number count angular power spectra

# Description

The code is a fully-operational modified version of the publicly available package code CosmoSIS ([Zunts et al. 2015](https://arxiv.org/pdf/1409.3409)). In this modular cosmological parameter estimation suite for galaxy number count Limber approximated angular power spectra, we allow for single and multiple tracers, including density fluctuations, redshift-space distortions, and weak lensing magnification.

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

the weak lensing magnification.
