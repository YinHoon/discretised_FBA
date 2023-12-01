"""
A branched pathway toy model for simulation
"""

METABOLITES = ['X[e]', 'X[c]', 'Y[c]', 'Energy', 'Toxin']
EXCHANGE_REACTION = {'EX_X': ' --> X[e]'}
TRANSPORT_REACTION = {'Transport_X': 'X[e] --> X[c]'}
INTRACELLULAR_REACTIONS = {
    'Synthesise_Y': 'X[c] --> Y[c]',
    'Primary': 'Y[c] --> 4 Energy',
    'Secondary': 'Y[c] --> 3 Energy + Toxin',
    'Demand_energy': 'Energy --> ',
    'Demand_toxin': 'Toxin --> '
    }
PROTEIN_REACTION_RULE = {
    'Transport_X': {'Protein1': 1.0},
    'Synthesise_Y': {'Protein2': 1.0},
    'Primary': {'Protein3': 1.0},
    'Secondary': {'Protein4': 2.0}
    }
CELLULAR_PROTEIN_CONCENTRATION = [
    ('Protein1', 100),
    ('Protein2', 100), 
    ('Protein3', 100),
    ('Protein4', 100)
]
CELL_SHAPES = [
    {'ShapeID': '6x6', 'Width': 6, 'Length': 6},
    {'ShapeID': '4x9', 'Width': 4, 'Length': 9},
    {'ShapeID': '3x12', 'Width': 3, 'Length': 12},
    {'ShapeID': '2x18', 'Width': 2, 'Length': 18},
    {'ShapeID': '1x36', 'Width': 1, 'Length': 36},
]