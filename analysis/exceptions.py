#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Custom exception classes."""

__all__ = ['NoCSPData', 'UnobservedFeature']


class NoCSPData(Exception):
    """There is no CSP published t0 or E(B - V) value for this target"""
    pass


class UnobservedFeature(Exception):
    """Referencing a spectral feature that is at least partially unobserved."""
    pass
