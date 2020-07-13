The Fisher sampler
--------------------------------------------------------------------

Fisher matrix calculation

+--------------+------------------------------------------+
| | Name       | | fisher                                 |
+--------------+------------------------------------------+
| | Version    | | 1.0                                    |
+--------------+------------------------------------------+
| | Author(s)  | | CosmoSIS Team                          |
+--------------+------------------------------------------+
| | URL        | | https://bitbucket.org/joezuntz/cosmosis|
+--------------+------------------------------------------+
| | Citation(s)|                                          |
+--------------+------------------------------------------+
| | Parallelism| | embarrassing                           |
+--------------+------------------------------------------+



The Fisher information matrix characterizes the curvature of a distribution, typically at its peak.  For multi-dimensional parameter spaces it can be given in  matrix form, with each element representing the curvature between two parameters.

Fisher matrices can be used to approximate the full likelihood shape one would derive from a sampling process, but much faster as only a handful of likelihood evaluations are needed.  This approximation is exact for perfectly Gaussian posteriors but not otherwise; in particular for distributions with "Banana-like" degeneracies or cut-offs near the peak it is a poor approximation.

For a general distribution the Fisher matrix can be rather fiddly to calculate;  this sampler does not do that.  Instead it assumes a Gaussian likelihood, in which case the fisher matrix can be computed from to the derivatives of the observable quantities :math:`v` with respect to the parameters :math:`p`, and their covariance matrix :math:`C` :

.. math::

    F_{ij} = \sum_{mn} \frac{\partial v_m}{\partial p_i} [C^{-1}]_{mn} \frac{\partial v_n}{\partial p_j}



There is an additional term that arises when the covariance matrix C depends on the parameters p, which we do not currently calculate here, as it is usually fixed in  cosmology. We plan to implement this shortly, however.

The CosmoSIS fisher sampler is calculated around the central value provided in the values file; no optimization is done before running.

Unlike most other CosmoSIS samplers which depend only likelihood of a parameter set, the fisher sampler requires the predicted observables from a pipeline too.  They are expected to be saved to a section of the data block, called "data_theory", with  keys, name+"_theory" and name+"_inverse_covariance" where "name" is the name of the likelihood, for example, to save a cmb data vector you could save "cmb_theory" and "cmb_inverse_covariance" in the "data_vector" section.

Gaussian likelihoods implemented in CosmoSIS save these sections automatically. Your own likelihoods can use the CosmoSIS gaussian likelihood superclass to do the same, or you can manually save name+"_theory" and name+"_inverse_covariance" for any data that  you add.



Installation
============

No special installation required; everything is packaged with CosmoSIS




Parameters
============

These parameters can be set in the sampler's section in the ini parameter file.  
If no default is specified then the parameter is required. A listing of "(empty)" means a blank string is the default.

+------------+--------+-----------------------------------------------------------+----------+
| | Parameter| | Type | | Meaning                                                 | | Default|
+------------+--------+-----------------------------------------------------------+----------+
| | step_size| | float| | The size, as a fraction of the total parameter range, of| | 0.01   |
|            |        | | steps to use in the derivative calculation. You should  |          |
|            |        | | investigate stability wrt this.                         |          |
+------------+--------+-----------------------------------------------------------+----------+
