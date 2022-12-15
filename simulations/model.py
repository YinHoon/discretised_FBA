"""
A toy model for simulation
"""

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