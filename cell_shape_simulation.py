""" Design of simulation experiments

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Date: 07-11-2022
:License: MIT

"""
import itertools
from discretised_fba import DiscretisedCell


METABOLITES = ['A[e]', 'A[c]', 'B[c]', 'Biomass']
EXCHANGE_REACTION = {'EX_A': ' --> A[e]'}
TRANSPORT_REACTION = {'Transport_A': 'A[e] --> A[c]'}
INTRACELLULAR_REACTIONS = {
    'Synthesise_B': 'A[c] --> B[c]',
    'Synthesise_biomass': 'B[c] --> Biomass',
    'Biomass_reaction': 'Biomass --> '
    }
PROTEIN_REACTION_RULE = {
    'Transport_A': {'Protein1': 1.0},
    'Synthesise_B': {'Protein2': 1.0},
    'Synthesise_biomass': {'Protein3': 1.0}
    }
CELLULAR_PROTEIN_CONCENTRATION = [
    ('Protein1', 100),
    ('Protein2', 100), 
    ('Protein3', 100),
]
CELL_SHAPES = [
    {'ShapeID': '6x6', 'Width': 6, 'Length': 6},
    {'ShapeID': '4x9', 'Width': 4, 'Length': 9},
    {'ShapeID': '3x12', 'Width': 3, 'Length': 12},
    {'ShapeID': '2x18', 'Width': 2, 'Length': 18},
    {'ShapeID': '1x36', 'Width': 1, 'Length': 36},
]
ENZYME_DISTRIBUTION = ['uniform', 'moderate_inside_out', 'steep_inside_out', 
    'moderate_outside_in','steep_outside_in', 'random']


for shape in CELL_SHAPES:
    cell = DiscretisedCell(shape['ShapeID'], shape['Width'], shape['Length'])
    cell.create_reactions(METABOLITES[0], METABOLITES[1:], EXCHANGE_REACTION,
        INTRACELLULAR_REACTIONS, 'Biomass_reaction')
    cell.create_transport_reactions(['A[c]'], TRANSPORT_REACTION)
    cell.create_diffusion(['A[c]', 'B[c]'])
    
    permutations = list(itertools.product(ENZYME_DISTRIBUTION, repeat=3))
    for i in permutations:
        # remove gradient distribution for transport protein 
        if i[0] not in ['inside_out', 'outside_in']:
            for ind, protein in enumerate(CELLULAR_PROTEIN_CONCENTRATION):
                distribution = i[ind]
                all_regions = False if ind == 0 else True # transport protein
                random_distribution = True if distribution == 'random' else False
                if distribution == 'moderate_inside_out':
                    gradient = -0.5
                elif distribution == 'steep_inside_out':
                    gradient = -1.0
                elif distribution == 'moderate_outside_in':
                    gradient = 0.5
                elif distribution == 'steep_outside_in':
                    gradient = 1.0
                else:
                    gradient = 0.0    
                cell.distribute_enzyme(protein[0], protein[1], gradient=gradient, 
                    random_distribution=random_distribution, all_regions=all_regions)
        synth_rxns, trans_rxns = cell.calc_enzymatic_bounds(PROTEIN_REACTION_RULE)
        rxn_bounds = {}
        for rxn in cell.model.reactions:
            if rxn.id in synth_rxns:
                rxn_bounds[rxn.id] = (0, synth_rxns[rxn.id])
            elif rxn.id in trans_rxns:
                rxn_bounds[rxn.id] = (0, trans_rxns[rxn.id])
            else:
                rxn_bounds[rxn.id] = (0, 1000.)
        cell.set_bounds(rxn_bounds)
        solution = cell.solve()            
        print(f'{shape["ShapeID"]}:{i}')
        print(solution.objective_value)
