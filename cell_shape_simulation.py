""" Design of simulation experiments

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Date: 07-11-2022
:License: MIT

"""
METABOLITES = ['A[e]', 'A[c]', 'B[c]', 'Biomass']

EXCHANGE_REACTION = {'EX_A': ' --> A[e]'}

TRANSPORT_REACTION = {'Transport_A': 'A[e] --> A[c]'}

INTRACELLULAR_REACTIONS = {
	'Synthesise_B': 'A[c] --> B[c]',
	'Synthesise_biomass': 'B[c] --> Biomass'
	}

PROTEIN_REACTION_RULE = {
    'Transport_A': 'Protein1',
    'Synthesise_B': 'Protein2',
    'Synthesise_biomass': 'Protein3'}	

OBJECTIVE_REACTION = {'Biomass_reaction': 'Biomass --> '}

CELL_SHAPES = [
    {'ShapeID': '6x6', 'Width': 6, 'Length': 6},
    {'ShapeID': '3x12', 'Width': 3, 'Length': 12},
    {'ShapeID': '2x18', 'Width': 2, 'Length': 18},
    {'ShapeID': '1x36', 'Width': 1, 'Length': 36},
]

ENZYME_DISTRIBUTION = ['uniform', 'inside_out', 'outside_in', 'random']    