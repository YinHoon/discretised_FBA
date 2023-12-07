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
import numpy as np
import pandas as pd
import seaborn as sns
import sys

# Find code directory relative to our directory
THIS_DIR = dirname(__file__)
CODE_DIR = abspath(join(THIS_DIR, '..', '..', 'discretised_fba'))
sys.path.append(CODE_DIR)
from discretised_fba import DiscretisedCell

POPULATION_SIZE = 100
ENZYME_DISTRIBUTION = ['uniform', 'steep_inside_out', 'steep_outside_in', 'random']
DIFFUSING_METABOLITES = [['X[c]', 'Y[c]'], ['X[c]'], ['Y[c]'], []] 

def simulation():
    """ Run simulation

    Returns:
        :obj:`dict`: simulation results
    """
    # initialise results dictionary
    permutations = list(itertools.product(ENZYME_DISTRIBUTION, repeat=4))

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
                INTRACELLULAR_REACTIONS, {'Demand_energy': 1.0, 'Demand_toxin': -2.0})
            cell.create_transport_reactions(['X[c]'], TRANSPORT_REACTION)
            cell.create_diffusion(diffusion)        
            
            for distrib_name in results:
                distribution_combination = distrib_name.split(',')
                all_solutions = []
                
                for individual in range(POPULATION_SIZE):                       
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
                    all_solutions.append(cell.solve())
                results[distrib_name][diff_name].append({
                    'ShapeID': shape['ShapeID'],
                    'Aspect ratio': cell.aspect_ratio,
                    'Perimeter-to-area ratio': cell.peri_to_area,
                    'Simulation': all_solutions,
                    })

    return results            

def plot_cv(results, plot_name):
    """Plot coefficient of variation

    Args:
        results (:obj:`dict`): simulation results
        plot_name (:obj:`str`): name of plot
    """
    distribution_dictionary = {
            'uniform': '0', 
            'steep_inside_out': '\u2013', 
            'steep_outside_in': '+',
            'random': 'r'}
    cols = results.keys()
    colours = ['crimson', 'seagreen', 'cornflowerblue', 'orange']
    styles = ['dashed', 'solid', 'dashed', 'solid']
    fig, axes = plt.subplots(nrows=9, ncols=12, figsize=(10, 6), sharex=True, 
        sharey=True, constrained_layout=True)
    axes = axes.ravel()  # array to 1D
    fig.suptitle(plot_name)
    fig.supylabel('Coefficient of variation')
    fig.supxlabel('Perimeter-to-area ratio')
    cv = lambda x: np.std(x, ddof=1) / np.mean(x) * 100

    for col, ax in zip(cols, axes):
        for diff, colour, style in zip(sorted(results[col].keys()), colours, styles):
            value = results[col][diff]
            diff_name = diff if diff else 'No diffusion'
            x = []
            y = []            
            for i in value:                
                all_cells_data = []                
                for cell in i['Simulation']:
                    total_primary_flux = sum([cell.fluxes[j] \
                        for j in cell.fluxes.keys() \
                        if 'Primary' in j])
                    total_secondary_flux = sum([cell.fluxes[j] \
                        for j in cell.fluxes.keys() \
                        if 'Secondary' in j])
                    if plot_name == 'Total energy produced':
                        all_cells_data.append(4*total_primary_flux + 3*total_secondary_flux)
                    elif plot_name == 'Objective value':
                        all_cells_data.append(cell.objective_value)
                x.append(i['Perimeter-to-area ratio'])
                y.append(cv(all_cells_data))                
            
            ax.plot(x, y, colour, label=diff_name, linestyle=style, linewidth=2, alpha=0.35)
        
        subplot_title = ''.join([distribution_dictionary[i] for i in col.split(',')])
        ax.set_title(subplot_title, size=8)
        ax.xaxis.set_tick_params(labelsize=8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.set_ylim([0, 15])    
        line, label = ax.get_legend_handles_labels()
    
    fig.legend(line, label, loc='lower right', fontsize=9, bbox_to_anchor=(1.0, 0.))
    #fig.tight_layout()
    for i in range(-7, 0, 1):
        axes.flat[i].set_visible(False)
    dest_figname = abspath(join(THIS_DIR, '..', '..', 'results', 'branch_pathway',
        f'{plot_name}_cv.png'))
    fig.savefig(dest_figname, dpi=1200)


if __name__ == '__main__':
    results = simulation()
    plot_cv(results, 'Total energy produced')
    plot_cv(results, 'Objective value')