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
        self.assertEqual(len(cell.regions), 2)
        self.assertEqual(all(isinstance(i, Region) for j in cell.regions for i in j), True)
        self.assertEqual(sorted([i.id for j in cell.regions for i in j]), 
            sorted(['0,0', '0,1', '0,2', '1,0', '1,1', '1,2']))

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
        obj_rxn_id = 'R_biomass_ex'

        cell = DiscretisedCell('good_cell', 2, 1)
        cell.create_reactions(extra_mets, intra_mets, extra_rxns, intra_rxns, 
            obj_rxn_id)

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
        self.assertEqual(sorted(list(cell.regions[0,0].synthesis_reactions.keys())), sorted(intra_rxns))
        self.assertEqual(sorted(list(cell.regions[0,1].synthesis_reactions.keys())), sorted(intra_rxns))

        self.assertEqual(str(cell.model.objective.expression), 
            '1.0*R_biomass_ex_0,0 - 1.0*R_biomass_ex_0,0_reverse_db2c7 + ' +
            '1.0*R_biomass_ex_0,1 - 1.0*R_biomass_ex_0,1_reverse_f5050')

        bad_cell1 = DiscretisedCell('bad_cell1', 1, 2)
        obj_rxn_id = 'bad_obj_reaction'
        with self.assertRaises(KeyError) as error:
            bad_cell1.create_reactions(extra_mets, intra_mets, extra_rxns, intra_rxns, 
                obj_rxn_id)
            self.assertEqual(error.exception.message, 
                'obj_rxn_id is not an existing reaction')        

    def test_create_transport_reactions(self):
        intra_mets = ['A[c]']
        extra_intra_rxns = {'R_trans': 'A[e] --> A[c]'}
        
        cell = DiscretisedCell('good_cell', 4, 3)
        cell.create_transport_reactions(intra_mets, extra_intra_rxns)

        test_reactions = {
            'R_trans_0,0': 'A[e] --> A[c]_0,0', 
            'R_trans_0,1': 'A[e] --> A[c]_0,1',
            'R_trans_0,2': 'A[e] --> A[c]_0,2',
            'R_trans_0,3': 'A[e] --> A[c]_0,3',
            'R_trans_1,0': 'A[e] --> A[c]_1,0',
            'R_trans_1,3': 'A[e] --> A[c]_1,3',
            'R_trans_2,0': 'A[e] --> A[c]_2,0',
            'R_trans_2,1': 'A[e] --> A[c]_2,1',
            'R_trans_2,2': 'A[e] --> A[c]_2,2',
            'R_trans_2,3': 'A[e] --> A[c]_2,3',
            }
        self.assertDictEqual({i.id: i.reaction for i in cell.model.reactions}, test_reactions)
        for i in [(0,0), (0,1), (0,2), (0,3), (1,0), (1,3), (2,0), (2,1), (2,2), (2,3)]:
            self.assertEqual(sorted(list(cell.regions[i].transport_reactions.keys())), sorted(extra_intra_rxns))

    def test_create_diffusion(self):
        
        intra_mets = ['A[c]']
        cell = DiscretisedCell('test_cell', 3, 4)
        cell.create_diffusion(intra_mets)

        diffusion_region_0_0 = {
            'A[c]_0,0_to_0,1': 'A[c]_0,0 --> A[c]_0,1',
            'A[c]_0,0_to_1,0': 'A[c]_0,0 --> A[c]_1,0',
            'A[c]_0,0_to_1,1': 'A[c]_0,0 --> A[c]_1,1',
            'A[c]_0,1_to_0,0': 'A[c]_0,1 --> A[c]_0,0',
            'A[c]_1,0_to_0,0': 'A[c]_1,0 --> A[c]_0,0',
            'A[c]_1,1_to_0,0': 'A[c]_1,1 --> A[c]_0,0',
            }
        self.assertEqual({k:v.reaction for k, v in cell.regions[0,0].diffusions.items()}, 
            diffusion_region_0_0)             

    def test_get_neighbouring_regions(self):
        cell = DiscretisedCell('good_cell', 3, 4)

        self.assertEqual(sorted(cell.get_neighbouring_regions(0,0)), 
            sorted([(0,1), (1,0), (1,1)]))
        self.assertEqual(sorted(cell.get_neighbouring_regions(0,1)), 
            sorted([(0,0), (0,2), (1,0), (1,1), (1,2)]))
        self.assertEqual(sorted(cell.get_neighbouring_regions(1,0)), 
            sorted([(0,0), (0,1), (1,1), (2,0), (2,1)]))
        self.assertEqual(sorted(cell.get_neighbouring_regions(1,1)), 
            sorted([(0,0), (0,1), (0,2), (1,0), (1,2), (2,0), (2,1), (2,2)]))

    def test_distribute_enzyme(self):
        # uniform distribution to all regions
        cell = DiscretisedCell('cell', 3, 4)
        cell.distribute_enzyme('enzyme1', 120)
        for i in cell.regions:
            for j in i:
                self.assertEqual(j.enzymes['enzyme1'], 10)
        self.assertAlmostEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 120, delta=1e-12)
        # high gradient distribution: decreasing outside in
        cell = DiscretisedCell('cell', 6, 6)
        cell.distribute_enzyme('enzyme1', 550, gradient=-1)
        first_layer = [(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(0,1),(5,1),(0,2),
            (5,2),(0,3),(5,3),(0,4),(5,4),(0,5),(1,5),(2,5),(3,5),(4,5),(5,5)]
        for i in first_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 18.75)
        second_layer = [(1,1),(2,1),(3,1),(4,1),(1,2),(4,2),(1,3),(4,3),
            (1,4),(2,4),(3,4),(4,4)]
        for i in second_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 12.5)
        third_layer = [(2,2),(2,3),(3,2),(3,3)]
        for i in third_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 6.25)    
        self.assertEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 550)
        # low gradient distribution: decreasing outside in
        cell = DiscretisedCell('cell', 6, 6)
        cell.distribute_enzyme('enzyme1', 1060, gradient=-0.1)
        first_layer = [(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(0,1),(5,1),(0,2),
            (5,2),(0,3),(5,3),(0,4),(5,4),(0,5),(1,5),(2,5),(3,5),(4,5),(5,5)]
        for i in first_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 30)
        second_layer = [(1,1),(2,1),(3,1),(4,1),(1,2),(4,2),(1,3),(4,3),
            (1,4),(2,4),(3,4),(4,4)]
        for i in second_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 29)
        third_layer = [(2,2),(2,3),(3,2),(3,3)]
        for i in third_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 28)    
        self.assertEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 1060)
        # high gradient distribution: increasing outside in
        cell = DiscretisedCell('cell', 6, 6)
        cell.distribute_enzyme('enzyme1', 560, gradient=1)
        first_layer = [(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(0,1),(5,1),(0,2),
            (5,2),(0,3),(5,3),(0,4),(5,4),(0,5),(1,5),(2,5),(3,5),(4,5),(5,5)]
        for i in first_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 10)
        second_layer = [(1,1),(2,1),(3,1),(4,1),(1,2),(4,2),(1,3),(4,3),
            (1,4),(2,4),(3,4),(4,4)]
        for i in second_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 20)
        third_layer = [(2,2),(2,3),(3,2),(3,3)]
        for i in third_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 30)    
        self.assertEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 560)
        # low gradient distribution: increasing outside in
        cell = DiscretisedCell('cell', 6, 6)
        cell.distribute_enzyme('enzyme1', 460, gradient=0.5)
        first_layer = [(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(0,1),(5,1),(0,2),
            (5,2),(0,3),(5,3),(0,4),(5,4),(0,5),(1,5),(2,5),(3,5),(4,5),(5,5)]
        for i in first_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 10)
        second_layer = [(1,1),(2,1),(3,1),(4,1),(1,2),(4,2),(1,3),(4,3),
            (1,4),(2,4),(3,4),(4,4)]
        for i in second_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 15)
        third_layer = [(2,2),(2,3),(3,2),(3,3)]
        for i in third_layer:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 20)    
        self.assertEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 460)
        # random distribution to all regions
        cell = DiscretisedCell('cell', 3, 4)
        cell.distribute_enzyme('enzyme1', 120, random_distribution=True)
        self.assertAlmostEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 120, delta=1e-12)
        # uniform distribution to outermost layer
        cell = DiscretisedCell('cell', 3, 4)
        cell.distribute_enzyme('enzyme1', 120, all_regions=False)
        outer_coordinate = [(0,0), (0,1), (0,2), (1,0), (1,2), (2,0), (2,2), (3,0), (3,1), (3,2)]
        for i in outer_coordinate:
            self.assertEqual(cell.regions[i].enzymes['enzyme1'], 12)
        self.assertAlmostEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 120, delta=1e-12)  
        # random distribution to outermost layer
        cell = DiscretisedCell('cell', 3, 4)
        cell.distribute_enzyme('enzyme1', 120, random_distribution=True, all_regions=False)
        self.assertAlmostEqual(sum([cell.regions[i].enzymes['enzyme1'] for i in outer_coordinate]), 120, delta=1e-12)  
        self.assertAlmostEqual(sum([j.enzymes['enzyme1'] for i in cell.regions for j in i]), 120, delta=1e-12)

    def test_calc_enzymatic_bounds(self):
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
        enzyme_rxn_assoc = {'R_syn': {'e1': 0.5,'e2': 1.0}, 'R_trans': {'e1': 0.2}}

        cell = DiscretisedCell('good_cell', 2, 1)
        cell.create_reactions(extra_mets, intra_mets, extra_rxns, intra_rxns, 
            obj_rxn_id)
        cell.create_transport_reactions(['A[c]'], extra_intra_rxns)
        cell.distribute_enzyme('e1', 10)
        cell.distribute_enzyme('e2', 100)
        synth_rxn_bounds, trans_rxn_bounds = cell.calc_enzymatic_bounds(enzyme_rxn_assoc)

        test_syn_rxn_bounds = {'R_syn_0,0': 52.5, 'R_syn_0,1': 52.5}
        test_trans_rxn_bounds = {'R_trans_0,0': 1.0, 'R_trans_0,1': 1.0}
        self.assertDictEqual(synth_rxn_bounds, test_syn_rxn_bounds)
        self.assertDictEqual(trans_rxn_bounds, test_trans_rxn_bounds)

    def test_set_bounds(self):
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

        cell = DiscretisedCell('good_cell', 4, 3)
        cell.create_reactions(extra_mets, intra_mets, extra_rxns, intra_rxns, 
            obj_rxn_id)
        cell.create_transport_reactions(['A[c]'], extra_intra_rxns)
        cell.create_diffusion(['A[c]'])
        rxn_bounds = {}
        for x in cell.regions:
            for region in x:
                for rxn in region.synthesis_reactions.values():
                    rxn_bounds[rxn.id] = (0, 10.)
                for rxn in region.transport_reactions.values():
                    rxn_bounds[rxn.id] = (0, 10.)    
                for diff in region.diffusions.values():
                    rxn_bounds[diff.id] = (0, 5.)
        rxn_bounds['R_A_ex'] = (0, 10.)               
        cell.set_bounds(rxn_bounds)

        self.assertEqual({i.id:i.bounds for i in cell.model.reactions for k in intra_rxns.keys() if k in i.id},
            {i.id:(0, 10.) for i in cell.model.reactions for k in intra_rxns.keys() if k in i.id})
        self.assertEqual({i.id:i.bounds for i in cell.model.reactions for k in extra_intra_rxns.keys() if k in i.id},
            {i.id:(0, 10.) for i in cell.model.reactions for k in extra_intra_rxns.keys() if k in i.id})
        
        for x in cell.regions:
            for region in x:
                for diff in region.diffusions.values():
                    self.assertEqual(cell.model.reactions.get_by_id(diff.id).bounds, (0, 5.))         

    def test_solve(self):
        
        cell = DiscretisedCell('good_cell', 4, 3)
        cell.create_reactions(['A[e]'], ['A[c]'], {'R_A_ex': ' --> A[e]'}, {'R_syn': 'A[c] --> '}, 'R_syn')
        cell.create_transport_reactions(['A[c]'], {'R_trans': 'A[e] --> A[c]'})
        rxn_bounds = {rxn.id: (0, 20.) for rxn in cell.model.reactions}
        cell.set_bounds(rxn_bounds)
        solution = cell.solve()
        
        self.assertEqual(solution.objective_value, 20.)

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
        region.add_metabolite('m', met)
        self.assertEqual(region.metabolites, {'m': met})

    def test_add_synthesis_reaction(self):
        rxn = cobra.Reaction('R')
        region = Region(1, 2)
        region.add_synthesis_reaction('r', rxn)
        self.assertEqual(region.synthesis_reactions, {'r': rxn})

    def test_add_transport_reaction(self):
        rxn = cobra.Reaction('R')
        region = Region(1, 2)
        region.add_transport_reaction('r', rxn)
        self.assertEqual(region.transport_reactions, {'r': rxn})     

    def test_add_diffusion(self):
        rxn = cobra.Reaction('R')
        region = Region(1, 2)
        region.add_diffusion(rxn)
        self.assertEqual(region.diffusions, {'R': rxn})  

    def test_add_enzyme_concentration(self):
        region = Region(1,2)
        region.add_enzyme_concentration('protein1', 20.5)
        self.assertEqual(region.enzymes, {'protein1': 20.5})               
