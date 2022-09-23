""" Test methods for running discretised FBA

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Date: 15-09-2022
:License: MIT

"""

from discretised_fba import DiscretisedCell, Region
import cobra
import unittest


class CellTestCase(unittest.TestCase):
    def test_init(self):

        cell = DiscretisedCell('good_cell', 3, 2)
        self.assertEqual(cell.id, 'good_cell')
        self.assertEqual(cell.width, 3)
        self.assertEqual(cell.length, 2)
        self.assertEqual(cell.regions.size, 6)
        self.assertEqual(len(cell.regions), 3)
        self.assertEqual(all(isinstance(i, Region) for j in cell.regions for i in j), True)
        self.assertEqual(sorted([i.id for j in cell.regions for i in j]), 
            sorted(['0,0', '0,1', '1,0', '1,1', '2,0', '2,1']))

        with self.assertRaises(TypeError) as error:
            bad_cell1 = DiscretisedCell('bad_cell1', 0.2, 2)
            self.assertEqual(error.exception.message, 'width must be an integer')
            bad_cell2 = DiscretisedCell('bad_cell2', 2, 3.2)
            self.assertEqual(error.exception.message, 'length must be an integer')                

    def test_create_reactions(self):
        extra_mets = ['A[e]']
        intra_mets = ['A[c]', 'B[c]', 'Biomass[c]']
        extra_rxns = {'R_A_ex': ' --> A[e]'}        
        intra_rxns = {
            'R_syn': 'A[c] --> B[c]', # synthesis reaction
            'R_biomass': '2 B[c] --> Biomass[c]', # biomass reaction
            'R_biomass_ex': 'Biomass[c] -->  ' # biomass exchange reaction
            }
        extra_intra_rxns = {'R_trans': 'A[e] --> A[c]'}     
        obj_rxn_id = 'R_biomass_ex'

        cell = DiscretisedCell('good_cell', 1, 2)
        cell.create_reactions(extra_mets, intra_mets, extra_rxns, intra_rxns, 
            extra_intra_rxns, obj_rxn_id)

        test_metabolites = ['A[e]', 'A[c]_0,0', 'A[c]_0,1', 'B[c]_0,0', 
            'B[c]_0,1', 'Biomass[c]_0,0', 'Biomass[c]_0,1']
        self.assertEqual(sorted([i.id for i in cell.model.metabolites]), sorted(
            test_metabolites))

        test_reactions = {
            'R_A_ex': ' --> A[e]',
            'R_syn_0,0': 'A[c]_0,0 --> B[c]_0,0', # synthesis reaction
            'R_biomass_0,0': '2.0 B[c]_0,0 --> Biomass[c]_0,0', # biomass reaction
            'R_biomass_ex_0,0': 'Biomass[c]_0,0 --> ', # biomass exchange reaction
            'R_syn_0,1': 'A[c]_0,1 --> B[c]_0,1', # synthesis reaction
            'R_biomass_0,1': '2.0 B[c]_0,1 --> Biomass[c]_0,1', # biomass reaction
            'R_biomass_ex_0,1': 'Biomass[c]_0,1 --> ' # biomass exchange reaction
            }
        self.assertDictEqual({i.id: i.reaction for i in cell.model.reactions}, test_reactions)

        bad_cell1 = DiscretisedCell('bad_cell1', 1, 2)
        obj_rxn_id = 'bad_obj_reaction'
        with self.assertRaises(KeyError) as error:
            bad_cell1.create_reactions(extra_mets, intra_mets, extra_rxns, intra_rxns, 
                extra_intra_rxns, obj_rxn_id)
            self.assertEqual(error.exception.message, 
                'obj_rxn_id is not an existing reaction')
        

class RegionTestCase(unittest.TestCase):
    def test_init(self):
        region = Region(1, 2)
        self.assertEqual(region.id, '1,2')
        self.assertEqual(region.row_number, 1)
        self.assertEqual(region.column_number, 2)

        with self.assertRaises(TypeError) as error:
            bad_region1 = Region(0.2, 2)
            self.assertEqual(error.exception.message, 'row_number must be an integer')
            bad_region2 = Region(2, 3.2)
            self.assertEqual(error.exception.message, 'column_number must be an integer')

    def test_add_metabolite(self):
        met = cobra.Metabolite('M')
        region = Region(1, 2)
        region.add_metabolite(met)
        self.assertEqual(region.metabolites, {'M': met})

    def test_add_reaction(self):
        rxn = cobra.Reaction('R')
        region = Region(1, 2)
        region.add_reaction(rxn)
        self.assertEqual(region.reactions, {'R': rxn})                   
