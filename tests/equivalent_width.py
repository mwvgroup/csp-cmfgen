from unittest import TestCase

import numpy as np
from analysis import equivalent_width


class EWCalculation(TestCase):

    @classmethod
    def setUpClass(cls):
        """Create a top-hat absorption feature with a depth of 1"""

        cls.feat_width = 1.5
        cls.feat_wave = np.arange(-1, cls.feat_width + 1, .01)
        feat_bounds = [
            cls.feat_wave < 0,
            (0 <= cls.feat_wave) & (cls.feat_wave <= cls.feat_width),
            cls.feat_width < cls.feat_wave
        ]

        cls.feat_flux = np.piecewise(cls.feat_wave, feat_bounds, [1, 0, 1])

    def test_calc_ew(self):
        """Test equivalent_width.calc_pew"""

        ew = equivalent_width.calc_pew(
            self.feat_wave, self.feat_flux, 0, self.feat_width, lambda x: 1)

        self.assertEqual(self.feat_width, ew)
