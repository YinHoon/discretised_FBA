""" Test methods for running discretised FBA

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Date: 15-09-2022
:License: MIT

"""

from discretised_fba import DiscretisedCell
import unittest


class CellTestCase(unittest.TestCase):
    def test_init(self):

        cell = DiscretisedCell('good_cell', 1, 2)
        self.assertEqual(cell.id, 'good_cell')
        self.assertEqual(cell.width, 1)
        self.assertEqual(cell.length, 2)

        with self.assertRaises(TypeError) as error:
            bad_cell1 = DiscretisedCell('bad_cell1', 0.2, 2)
            self.assertEqual(error.exception.message, 'width must be an integer')
            bad_cell2 = DiscretisedCell('bad_cell2', 2, 3.2)
            self.assertEqual(error.exception.message, 'length must be an integer')

    	
