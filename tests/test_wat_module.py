from __future__ import print_function

import unittest
import warnings
import os

import numpy as np
np.set_printoptions(threshold=np.nan)

from numpy import array, float32
from numpy.testing import assert_almost_equal

from weppout import Wat

class Test_read(unittest.TestCase):
    def test_1ofe(self):
		import sys
		sys.path.append('../')
		
		wat_obj = Wat(r'data/1_wat_1ofe.txt')
		
		# Compare to the Dp data
#		with open(r'validation/test_1ofe_Dp.txt','w') as f:
#			f.write(repr(wat_obj.data['Dp (mm)']))
		validation_dp = eval(open(r'validation/test_1ofe_Dp.txt').read())
		assert_almost_equal(wat_obj.data['Dp (mm)'], validation_dp)
		
		self.assertEqual(wat_obj.num_ofes, 1)
		self.assertEqual(wat_obj.data['Dp (mm)'].shape, (365,1))
		self.assertEqual(wat_obj.data['J'].shape, (365,1))
		self.assertEqual(wat_obj.data['Area (m^2)'].shape, (1,1))
		self.assertEqual(wat_obj.num_ofes, 1)
		
    def test_3ofe(self):
		import sys
		sys.path.append('../')

		wat_obj = Wat(r'data/1_wat_3ofe.txt')
		
		# Compare to the Dp data
#		with open(r'validation/test_3ofe_Dp.txt','w') as f:
#			f.write(np.array_repr(wat_obj.data['Dp (mm)']))
		validation_dp = eval(open(r'validation/test_3ofe_Dp.txt').read())
		assert_almost_equal(wat_obj.data['Dp (mm)'], validation_dp)
		
		self.assertEqual(wat_obj.num_ofes, 3)
		self.assertEqual(wat_obj.data['Dp (mm)'].shape, (365,3))
		self.assertEqual(wat_obj.data['J'].shape, (365,1))
		self.assertEqual(wat_obj.data['Area (m^2)'].shape, (1,3))
		
			
def suite():
    return unittest.TestSuite((
            unittest.makeSuite(Test_read)
                              ))

if __name__ == "__main__":
    # run tests
    runner = unittest.TextTestRunner()
    runner.run(suite())
