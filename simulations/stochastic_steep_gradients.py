""" Simulations where at least one enzyme is randomly distributed

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Date: 07-11-2022
:License: MIT

"""

from model import *
from openpyxl import Workbook
from os.path import dirname, abspath, join
import collections
import itertools
import matplotlib.pyplot as plt
import random
import sys

# Find code directory relative to our directory
THIS_DIR = dirname(__file__)
CODE_DIR = abspath(join(THIS_DIR, '..', 'discretised_fba'))
sys.path.append(CODE_DIR)
from discretised_fba import DiscretisedCell

POPULATION_SIZE = 100
ENZYME_DISTRIBUTION = ['uniform', 'steep_inside_out', 'steep_outside_in', 'random']
DIFFUSING_METABOLITES = [['A[c]', 'B[c]'], ['A[c]'], ['B[c]'], []]  

# initialise results dictionary
permutations = list(itertools.product(ENZYME_DISTRIBUTION, repeat=3))

# remove gradient distribution for transport protein        
results = {','.join(i): {' and '.join(j): [] for j in DIFFUSING_METABOLITES} \
    for i in permutations if i[0] not in ['steep_inside_out', 'steep_outside_in'] \
    and 'random' in i}

# Run simulations
for diffusion in DIFFUSING_METABOLITES:
    diff_name = ' and '.join(diffusion)
    for shape in CELL_SHAPES:
        cell = DiscretisedCell(shape['ShapeID'], shape['Width'], shape['Length'])
        cell.create_reactions(METABOLITES[0], METABOLITES[1:], EXCHANGE_REACTION,
            INTRACELLULAR_REACTIONS, 'Biomass_reaction')
        cell.create_transport_reactions(['A[c]'], TRANSPORT_REACTION)
        cell.create_diffusion(diffusion)        
        
        for distrib_name in results:
            distribution_combination = distrib_name.split(',')
            all_solutions = []
            random.seed(1)
            for individual in range(POPULATION_SIZE):
                seed_value = random.randint(1,100)       
                for ind, protein in enumerate(CELLULAR_PROTEIN_CONCENTRATION):
                    distribution = distribution_combination[ind]
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
                        random_distribution=random_distribution, random_seed=seed_value, 
                        all_regions=all_regions)
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
                all_solutions.append(cell.solve())
            results[distrib_name][diff_name].append({
                'ShapeID': shape['ShapeID'],
                'Aspect ratio': cell.aspect_ratio,
                'Perimeter-to-area ratio': cell.peri_to_area,
                'Simulation': all_solutions,
                })
            
# Write results to an excel file
dest_filename = abspath(join(THIS_DIR, '..', 'results', 'stochastic_steep_gradients.xlsx'))
wb = Workbook()
sheets = {}

column_index = {distrib_name: ind+1  for ind, distrib_name in enumerate(results)}
for diffusion in DIFFUSING_METABOLITES:
    diff_name = ' and '.join([i[0:i.index('[')] for i in diffusion])
    diff_name = diff_name if diff_name else 'No diffusion'
    for shape in CELL_SHAPES:
        sheet_name = f"{shape['ShapeID']}({diff_name})"
        ws = wb.create_sheet(title=sheet_name)    
        sheets[sheet_name] = ws
        for distrib_name, ind in column_index.items():
            _ = ws.cell(column=ind, row=1, 
                        value=distrib_name)

for distrib_name, data in results.items():
    col = column_index[distrib_name]
    for diff, details in data.items():
        diffusion_list = diff.split(' and ')
        diff_name = ' and '.join([i[0:i.index('[')] for i in diffusion_list if i])
        diff_name = diff_name if diff_name else 'No diffusion'
        for shape in details:
            sheet_name = f"{shape['ShapeID']}({diff_name})"
            ws = sheets[sheet_name]
            row = 2
            for individual in shape['Simulation']:
                _ = ws.cell(column=col, row=row, 
                        value=individual.objective_value)
                row += 1        
        
wb.save(filename = dest_filename)
