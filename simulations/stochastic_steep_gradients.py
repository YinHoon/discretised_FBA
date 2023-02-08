""" Simulations where at least one enzyme is randomly distributed

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Date: 07-11-2022
:License: MIT

"""

from model import *
from openpyxl import Workbook
from os.path import dirname, abspath, join
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
import collections
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
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

# Plot kernel density
distribution_dictionary = {
        'uniform': '0', 
        'steep_inside_out': '\u2013', 
        'steep_outside_in': '+',
        'random': 'r'}
cols = results.keys()
colours = ['r', 'g', 'b', 'k']

for shape in CELL_SHAPES:
    fig, axes = plt.subplots(nrows=4, ncols=6, figsize=(10, 6), sharey=True)
    axes = axes.ravel()  # array to 1D
    fig.suptitle(f'Cell geometry: {shape["ShapeID"]}')
    for col, ax in zip(cols, axes):
        for diff, colour in zip(sorted(results[col].keys()), colours):
            value = results[col][diff]
            diff_name = diff if diff else 'No diffusion'
            data = [j.objective_value for i in value for j in i['Simulation'] \
                if i['ShapeID']==shape['ShapeID']]
            if np.var(data) < 1e-05:
                ax.axvline(np.mean(data), color=colour, label=diff_name)
            else:
                sns.set_style('white')
                sns.kdeplot(data, ax=ax, color=colour, label=diff_name)
        subplot_title = ''.join([distribution_dictionary[i] for i in col.split(',')])
        ax.title.set_text(subplot_title)
        ax.set_ylim([0, 0.3])
        ax.set_xlim([0, 120])
        line, label = ax.get_legend_handles_labels()
    fig.legend(line, label, loc='lower right', bbox_to_anchor=(1.0, 0.))
    fig.tight_layout()
    axes.flat[-1].set_visible(False)
    dest_figname = abspath(join(THIS_DIR, '..', 'results', f'{shape["ShapeID"]}.tiff'))
    fig.savefig(dest_figname, dpi=1200)

# Perform multivariate regression analysis on population data
column_headings = ['Protein1', 'Protein2', 'Protein3', 
    'Diffusion', 'Perimeter_to_area_ratio', 'Mean_growth', 
    'Standard_deviation_growth']
numerical_X_variables = ['Perimeter_to_area_ratio']
categorical_X_variables = ['Protein1', 'Protein2', 'Protein3', 'Diffusion']
Y_variables = ['Mean_growth', 'Standard_deviation_growth']
flatten_results = []
for distrib_name, v1 in results.items():
    distrib = distrib_name.split(',')
    for diffusion_type, v2 in v1.items():
        dif = diffusion_type if diffusion_type else 'No diffusion'
        for shape_sim in v2:
            data = [j.objective_value for j in shape_sim['Simulation']]
            if np.var(data) < 1e-05:
                data_entry = distrib + [dif]
                data_entry.append(shape_sim['Perimeter-to-area ratio'])            
                data_entry.append(np.mean(data))  
                data_entry.append(np.std(data))
                flatten_results.append(data_entry)        

df = pd.DataFrame(flatten_results, columns=column_headings)
# Standardize numerical variables
X_scaled = StandardScaler().fit_transform(df[numerical_X_variables])
X_scaled_df = pd.DataFrame(X_scaled, columns=numerical_X_variables)
Y_scaled = StandardScaler().fit_transform(df[Y_variables])
Y_scaled_df = pd.DataFrame(Y_scaled, columns=Y_variables)

# Convert categorical variables to dummy variables
categorical_df = pd.get_dummies(df[categorical_X_variables])
# Prepare dataframe for X and Y
X = pd.concat([X_scaled_df, categorical_df], axis=1)
Y_growth = Y_scaled_df[['Mean_growth']]
Y_std = Y_scaled_df[['Standard_deviation_growth']]

# Regression without interaction terms
model_growth = sm.OLS(Y_growth, sm.add_constant(X)).fit()
model_std = sm.OLS(Y_std, sm.add_constant(X)).fit()

# Regression with 2nd order interaction terms
interaction = PolynomialFeatures(degree=[2,2], include_bias=False, interaction_only=False)
transform_X = interaction.fit_transform(X)
interaction_terms = interaction.get_feature_names_out()
X_interaction2 = pd.DataFrame(transform_X, columns=interaction_terms)
model_growth_inter2 = sm.OLS(Y_growth, sm.add_constant(X_interaction2)).fit()
model_std_inter2 = sm.OLS(Y_std, sm.add_constant(X_interaction2)).fit()

# Regression with 3rd order interaction terms
interaction = PolynomialFeatures(degree=[3,3], include_bias=False, interaction_only=False)
transform_X = interaction.fit_transform(X)
interaction_terms = interaction.get_feature_names_out()
X_interaction3 = pd.DataFrame(transform_X, columns=interaction_terms)
model_growth_inter3 = sm.OLS(Y_growth, sm.add_constant(X_interaction3)).fit()
model_std_inter3 = sm.OLS(Y_std, sm.add_constant(X_interaction3)).fit()

dest_filename = abspath(join(THIS_DIR, '..', 'results', f'Multivariate_regression.csv'))
with open(dest_filename, 'w') as f:
    f.write(model_growth.summary().as_csv())
    f.write('\n\n')
    f.write(model_std.summary().as_csv())
    f.write('\n\n')
    f.write(model_growth_inter2.summary().as_csv())
    f.write('\n\n')
    f.write(model_std_inter2.summary().as_csv())
    f.write('\n\n')
    f.write(model_growth_inter3.summary().as_csv())
    f.write('\n\n')
    f.write(model_std_inter3.summary().as_csv())