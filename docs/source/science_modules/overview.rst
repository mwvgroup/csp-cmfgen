Section Outline
===============

This section outlines modules of the ``analysis`` package that are responsible for calculating and tabulating
scientific results. There is a dedicated module for each of the analysis steps outlined in the `Project Outline`_.

+-------------------------------+------------------------------------------------------------------------------------+
| Module Name                   | Description                                                                        |
+===============================+====================================================================================+
| ``analysis.spectra_chisq``    | Calculates the :math:`\chi^2` for observed spectra and a fiducial model SED.       |
+-------------------------------+------------------------------------------------------------------------------------+
| ``analysis.equivalent_width`` | Calculates the equivalent width of spectral features defined in `Folatelli+ 2004`_.|
+-------------------------------+------------------------------------------------------------------------------------+
| ``analysis.band_fitting``     | Independently fits separate filters of observed light-curves.                      |
+-------------------------------+------------------------------------------------------------------------------------+
| ``analysis.lc_colors``        | Interpolates the color of each supernova as a function of phase.                   |
+-------------------------------+------------------------------------------------------------------------------------+

.. _Folatelli+ 2004: https://ui.adsabs.harvard.edu/abs/2004NewAR..48..623F/abstract
