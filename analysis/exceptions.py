#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``exceptions`` module defines custom exception classes specific to the processing of CSP data.

Usage Example
-------------

.. code-block:: python
   :linenos:

   from analysis.exceptions import NoCSPData, UnobservedFeature
   
   raise NoCSPData('No data was available for a given target.')

   raise UnobservedFeature('A given feature is not spanned by a given spectrum.')

Function Documentation
----------------------
"""

__all__ = ['NoCSPData', 'UnobservedFeature']


class NoCSPData(Exception):
    """There is no CSP published t0 or E(B - V) value for this target"""
    pass


class UnobservedFeature(Exception):
    """Referencing a spectral feature that is at least partially unobserved."""
    pass
