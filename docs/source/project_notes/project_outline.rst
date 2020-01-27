Project Outline
===============

Overview
--------

Following the discovery of the accelerating expansion of the Universe (`Riess+ 1998`_, `Perlmutter+ 1999`_), Type Ia
supernovae (SNe Ia) have become a common tool for measuring the cosmological properties of the Universe
(`Betoule+ 2014`_, `DES 2019`_). However, the underlying physics that drive SN explosions is still not well understood.
Although it is thought that SNe Ia result from the thermonuclear explosions of white dwarfs (`Hoyle+ 1960`_,
`Woosley+ 1986`_), ongoing SN surveys have observed an increasing diversity in the SNe Ia population that is not well
described by existing models (see `Taubenberger+ 2017`_). The upcoming Legacy Survey of Space and Time (LSST) presents
an exciting new opportunity to answer long-standing questions by achieving the unprecedented: 20 TB of new data every
night, an entire survey of the southern night sky every three days, and on order a million new SNe discovered over a
ten years. However, with hundreds of new SNe discovered every night, it will be impossible to perform spectroscopic
follow-up observations for more than a small subset of discovered SNe.

Spectroscopic observations play a critical role in constraining SN models, primarily because their high wavelength
resolution promotes the identification of key physical processes that drive SN evolution. The impending shortage of
spectroscopic observations thus results in serious concerns over how well we will be able to constrain the fundamental
physics of SNe with LSST. This project is intended to address these concerns by 1) establishing the extent to which
existing SN models can be constrained using purely photometric observations, and 2) defining a prescription for
combining photometrically dominated data sets with limited spectroscopic follow-up to maximize constraints on the
parameter space of SN models.

.. _Riess+ 1998: https://ui.adsabs.harvard.edu/abs/1998AJ....116.1009R/abstract
.. _Perlmutter+ 1999: https://ui.adsabs.harvard.edu/abs/1999ApJ...517..565P/abstract
.. _Betoule+ 2014: https://ui.adsabs.harvard.edu/abs/2014A%26A...568A..22B/abstract
.. _DES 2019: https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.5329L/abstract
.. _Hoyle+ 1960: https://ui.adsabs.harvard.edu/abs/1960ApJ...132..565H/abstract
.. _Woosley+ 1986: https://ui.adsabs.harvard.edu/abs/1986BAAS...18.1016W/abstract
.. _Taubenberger+ 2017: https://ui.adsabs.harvard.edu/abs/2017hsn..book..317T/abstract

Key Steps
---------

As a fiducial set of SN models, this project will consider the four models outlined in (`Wilk+ 2017`_). These models
represent a range of SN masses, including two low mass models :math:`(\text{M} = 1.02 \text{ M}_\odot` and
:math:`1.04 \text{ M}_\odot)` an intermediate-mass model :math:`(\text{M} = 1.40 \text{ M}_\odot)`, and a super-massive
model :math:`(\text{M} = 1.70 \text{ M}_\odot)`. The considered models also allow us to consider a range of clumping)
factors in the explosion :math:`(f = 0.33, 0.25, 0.10)` in addition to the unclumped case :math:`(f=1`). This range in
simulated behavior ensures we will be able to consider cases where the models both perform well and struggle to
recreate the observations. To constrain each model, we will use spectroscopic and photometric observations of 134 SNe
taken by the Carnegie Supernova Project (CSP; `Krisciunas+ 2017`_). This comparison will be broken down into four key
stages. In the first two stages, we will use observations from the spectroscopic component of CSP to establish a
baseline for each model's performance. The latter two stages will then use CSP photometry to constrain the parameter
space of each model.

.. _Wilk+ 2017: https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.3187W/abstract
.. _Krisciunas+ 2017: https://ui.adsabs.harvard.edu/abs/2017AJ....154..211K/abstract

Stage 1
^^^^^^^

The first step in this project is to perform a direct comparison of the spectral energy distributions (SED) of each
model against the flux values of the observed spectra. In the ideal case, each model would be fitted directly to each
of the observed spectra. However, the spectral fitting process requires significant human oversight and is highly time
intensive. This makes fitting the entire set of CSP observed spectra (> 600 spectra) impractical. As an alternative,
we simulate a fiducial SED from each model at the phase of each observed spectrum. To limit the impact of background
noise, the :math:`\chi^2` of each model is then determined using a sliding window across the observed wavelength range.
The resulting values provide a quantitative description of how well each model represents CSP observations as a function
of phase (i.e., the time of the explosion) and wavelength. We can thus identify key phase and wavelength ranges where
the models diverge significantly from the observations and how the :math:`\chi^2` in these areas scales with model
parameters.

Stage 2
^^^^^^^

Having established an initial baseline for the performance of each model, we supplement it with a more detailed analysis
of key spectroscopic features. While an overall comparison of the SED is sufficient for identifying broad trends,
additional analysis is needed to determine detailed variations in critical chemical abundances of each SN. Specifically,
this comparison will be performed using 8 spectral features commonly used in the literature and outlined in
`Folatelli+ 2004`_. Not only does this allow for a comparison of our results with existing values in the literature, but
these also features represent chemical abundances that drive the evolutionary behavior of the observed spectra. By
identifying collections of modeled features that disagree with the observations, it is possible to isolate key
differences in the physical processes of each model that result in significant deviations from the observed behavior.

.. _Folatelli+ 2004: https://ui.adsabs.harvard.edu/abs/2004NewAR..48..623F/abstract

Stage 3
^^^^^^^

The third stage of this project shifts focus by using the models to fit the photometric observations directly.
Photometric observations are performed by using filters to divide observed light into broad-band wavelength bins.
When fitting models to photometric data, it is most common to fit all of the available data for a given SN
simultaneously. However, this project will take an alternative approach and separately fit the models to observations
taken in each bandpass. In doing so, model fits are constrained only by the way the photometric observations evolve in
that wavelength bin over time. This approach allows us to isolate how each model deviates from the observed behavior
in a given wavelength range and as a function of time, thus providing an analog of Stage 1 but with photometric data.

Stage 4
^^^^^^^

Although independently fitting each bandpass does provide valuable insight into each model’s behavior, it does not
account for how the ratio of flux between bandpasses, also known as the color, changes over time. This information
is useful because it provides insights into underlying physical behavior. For example, redder SNe, which have more
observed light at longer wavelengths than shorter ones, tend to have cooler temperatures and slower explosion
velocities. To incorporate this information, we will directly compare the color evolution predicted by each model
against the color values interpolated from the observed data. This tightens the existing constraints from Stage 3
while simultaneously providing a physical intuition for any observed deviations in the models.
