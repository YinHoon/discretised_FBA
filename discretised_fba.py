""" Methods for running discretised FBA

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Author: Alexandra Anamaria-Sorinca
:Date: 13-09-2022
:License: MIT

"""

import cobra
import numpy as np


class DiscretisedCell(object):
    """ A class of discretised cell objects. Each cell object has a 2D regular shape 
        defined by an MxN matrix that represents the cell's width (M) and length (N) 
    """
    def __init__(self, cell_id, width, length):
        """ 
        Args:
            id (:obj:`str`): cell ID
            width (:obj:`int`): cell width
            length (:obj:`int`): cell length

        Raises:
            :obj:`TypeError`: if the width and/or length is not an integer    
        """

        if not isinstance(width, int):
            raise TypeError('width must be an integer')
        if not isinstance(length, int):
            raise TypeError('length must be an integer')    

        self.id = cell_id
        self.width = width
        self.length = length
        self.regions = np.empty((self.width, self.length), dtype=np.object)
        self.model = cobra.Model()
        for i in range(self.width):
            for j in range(self.length):
                new_region = Region(i, j)
                self.regions[i, j] = new_region    

    def create_reactions(self, extra_mets, intra_mets, extra_rxns, intra_rxns, 
            obj_rxn_id):
        """ Create a Cobra model and add reactions to each region

        Args:
            extra_mets (:obj:`list` of :obj:`str`): list of IDs of metabolites 
                in the extracellular space
            intra_mets (:obj:`list` of :obj:`str`): list of IDs of metabolites 
                in the intracellular space
            extra_rxns (:obj:`dict`): dictionary of reaction ID as key and reaction
                string (e.g. '2 A --> B') as value, for reactions in the extracellular
                space
            intra_rxns (:obj:`dict`): dictionary of reaction ID as key and reaction
                string (e.g. '2 A --> B') as value, for reactions in the intracellular
                space
            obj_rxn_id (:obj:`str`): ID of objective reaction
            
        Raises:
            :obj:`KeyError`: if objective reaction ID is not in the
                extracellular, intracellular or transport reaction IDs            
        """
        if obj_rxn_id not in list(intra_rxns.keys()):
            raise KeyError(
                'obj_rxn_id is not an existing reaction')
        
        # create metabolites in the extracellular compartment
        metabolite_objs = []
        for i in extra_mets:
            metabolite_objs.append(cobra.Metabolite(i))
        self.model.add_metabolites(metabolite_objs)           
        
        # create reactions in the extracellular compartment
        reaction_objs = []
        for k in extra_rxns.keys():
            rxn_obj = cobra.Reaction(k)
            reaction_objs.append(rxn_obj)
        self.model.add_reactions(reaction_objs)    
        for rxn_obj in reaction_objs:
            rxn_obj.reaction = extra_rxns[rxn_obj.id]                    
        
        # create intracellular metabolites for each discrete region
        metabolite_objs = []
        for row in self.regions:
            for region in row:
                for i in intra_mets:
                    region_metabolite = cobra.Metabolite(f'{i}_{region.id}')
                    metabolite_objs.append(region_metabolite)
                    region.add_metabolite(region_metabolite)    
        self.model.add_metabolites(metabolite_objs)

        # create intracellular reactions for each discrete region
        for row in self.regions:
            for region in row:
                for k, v in intra_rxns.items():
                    region_rxn = cobra.Reaction(f'{k}_{region.id}')                    
                    region.add_reaction(region_rxn)
                    self.model.add_reactions([region_rxn])
                    new_rxn_str = v 
                    for i in intra_mets:
                        if i in new_rxn_str:
                            new_rxn_str = new_rxn_str.replace(i, f'{i}_{region.id}')
                    region_rxn.reaction = new_rxn_str

        # Set objective function
        objective_reactions = [rxn.flux_expression for rxn in self.model.reactions 
            if obj_rxn_id in rxn.id]
        objective_function = self.model.problem.Objective(sum(objective_reactions),
            direction='max')
        self.model.objective = objective_function

    def create_transport_reactions(self, intra_mets, extra_intra_rxns):
        """ Add nutrient transport reactions

        Args:
            intra_mets (:obj:`list` of :obj:`str`): list of IDs of intracellular 
                metabolites that are transported from the extracellular space
            extra_intra_rxns (:obj:`dict`): dictionary of reaction ID as key and 
                reaction string (e.g. '2 A --> B') as value, for transport reactions 
                between the extracellular and intracellular space   
        """
        # transport reactions between the extracellular and intracellular space
        x,y = np.ogrid[0:self.width, 0:self.length]
        # retrieve the inner layer
        inner_layer = (x>0)&(x<self.width-1)&(y>0)&(y<self.length-1)
        # create transport reaction to outer layer
        for region in self.regions[~inner_layer]:
            for k, v in extra_intra_rxns.items():
                region_rxn = cobra.Reaction(f'{k}_{region.id}')                    
                region.add_reaction(region_rxn)
                self.model.add_reactions([region_rxn])
                new_rxn_str = v 
                for i in intra_mets:
                    if i in new_rxn_str:
                        new_rxn_str = new_rxn_str.replace(i, f'{i}_{region.id}')
                region_rxn.reaction = new_rxn_str

            

class Region(object):
    """ A class of objects where each object represents a region on a cell.
        Each region occupies a space on a square lattice grid that is represented
        as a matrix. The position of each region is represented by its row and column 
        number on the matrix.
    """
    def __init__(self, row_number, column_number):
        """ 
        Args:
            row_number (:obj:`int`): row number
            column_number (:obj:`int`): column number

        Raises:
            :obj:`TypeError`: if the row and/or column number is not an integer    
        """
        if not isinstance(row_number, int):
            raise TypeError('row_number must be an integer')
        if not isinstance(column_number, int):
            raise TypeError('column_number must be an integer')    

        self.row_number = row_number
        self.column_number = column_number
        self.id = f'{row_number},{column_number}'
        self.metabolites = {}
        self.reactions = {}

    def add_metabolite(self, metabolite):
        """ Add a metabolite to the region

        Args:
            metabolite (:obj:`cobra.Metabolite`): a metabolite object
        """
        self.metabolites[metabolite.id] = metabolite

    def add_reaction(self, reaction):
        """ Add a metabolite to the region

        Args:
            metabolite (:obj:`cobra.Reaction`): a reaction object
        """
        self.reactions[reaction.id] = reaction    

#create the shape of the cell in a matrix form
#the row size * the column size = the number of cell compartiments
def run_model(r_size, c_size, enzyme):
    r_size = r_size
    c_size = c_size

    '''
    #if ij=1 then that part of a cell exists, otherwise ij=0
    #if ij=0 set the fluxes to be 0
    
    d = np.ones((int(r_size),int(c_size)))
    d[0,1] = '0'
    print(d)
    '''

    #define the reactions and their fluxes
    #set fluxes upper/lower bounds
    #metabolites create within the reaction models


    test_model = cobra.Model("test_model")

    # create empty dictionary for metabolite
    met = {}
    # create empty dictionary for reactions
    react = {}
    # create Ae and Re separately and add them to the met and react dict
    # Ae and Re are separated because they are from the external environment, which for simplicity is set to be neutral

    Ae = cobra.Metabolite('Ae')
    test_model.add_metabolites([Ae])
    met[Ae.id]=Ae
    rE = cobra.Reaction('rE')
    test_model.add_reactions([rE])
    rE.reaction = " --> Ae"
    react[rE.id]=rE

    #create names for the rest of the reactions, which are their id's
    #the react object in the dictionary is the reaction itself, defined by the metabolites
    #cobrapy has a tool which enables us to define the metabolites from the reaction....
    #instead of defining the metabolites manually

    for i in range(r_size):
        for j in range(c_size):
            #create name for metA ,B,C for compartment ij
            mets = ['A[c]','B[c]','Biomass[c]']
            for m in mets:
                met_id = m + str (i+1) + str (j+1)
                met_obj = cobra.Metabolite(met_id)
                test_model.add_metabolites([met_obj])
                met[met_id] = met_obj
            r_name = {'rAeAc': "Ae --> A[c]",
                      'rAcBc': "A[c] --> B[c]",
                      'rBcBio': "2 B[c] --> Biomass[c]",
                      'rBioE': "Biomass[c] -->  "}
            for k,v in r_name.items():
                react_id = k + str (i+1) + str(j+1)
                react_obj = cobra.Reaction(react_id)
                test_model.add_reactions([react_obj])
                react[react_id] = react_obj
                result = re.sub("]", "]" + str(i+1) + str(j+1), v)
                react_obj.reaction = result
                print(result)
    print(react)

    #enzymes can be distributed non-uniformly across the cell....
    #....which causes the geometric non-uniformity of the cell
    # in the excel workbook (appendix), there are a few ways to distribute the enzymes
    # this is currently done manually, each value is being inserted in the terminal after running the code

    enzyme_concentration = {'k_AeAc': 1000, 'k_AcBc':2000}
    distrib = {}
    for m in enzyme_concentration:
        for i in range(r_size):
            for j in range(c_size):
                distrib_id = m + str(i+1) + str(j+1)
                distrib[distrib_id]=enzyme[distrib_id]
        if not np.testing.assert_almost_equal(sum([v for k,v in distrib.items() if m in k]), 1, decimal=5, err_msg='"The concentration is not balanced"', verbose=True):
            pass
            print(sum([v for k,v in distrib.items() if m in k]))
    print(enzyme_concentration)
    print(distrib)

    # for each of the reactions within the dictionary, we need to set bounds for their fluxes
    #the lower bound will be 0, while the upper bound is determined by the concentration of the enzyme....
    #...which facilitates the reaction within a specific cell compartiment


    n = min(r_size,c_size)
    if n % 2 ==0:
        condition = n//2
    else:
        condition = n//2 + 1

    enzyme_concentration = {'k_AeAc': 1000, 'k_AcBc':2000}

    diffusion = {}

    for m, n in enzyme_concentration.items():
            for i in range(1, r_size+1):
                for j in range(1, c_size+1):
                    diffusion_id = m + str(i) + str(j)
                    for k, v in enzyme.items():
                        if diffusion_id == k:
                            if i == 1 or j == 1 or i == r_size or j == c_size:
                                diffusion[diffusion_id] = v * n
                            else:
                                for x in range(1, condition):
                                    if i != x and i != r_size - x + 1 and j != x and j != c_size - x + 1:
                                        diffusion[diffusion_id] = v * n * ((0.8) ** (x))

    print(diffusion)

    #some reactions have 2 indices and others have 3 (when the row numer has 2 digits)....
    #hence reactions id's can either have 7 or 8 characters
    #we are interested in comparing the last charachters from the reaction id with the last characters in enzyme id
    for k,v in react.items():
        for i, j in diffusion.items():
            if len(k) == 7:
                if k[-6:] == i[-6:]:
                    react[k].bounds=(0.0, j)
            else:
                if k[-7:] == i[-7:]:
                    react[k].bounds=(0.0, j)


    # let r/c be the row/column index of the cell part being deleted
    # if both the concentrations of the enzymes are being 0, then automatically the immediate neighbouring cell parts become external
    # notice that the deleted cell parts can only be external
    #this is just a suggestion for future work
    '''
    r = int(input("Deletion simulation row index "  + ":"))
    c = int(input("Deletion simulation column index "  + ":"))
    if r==0 or c==0:
        pass
    elif r!=1 and r!=r_size and c!=1 and c!=c_size:
        raise Exception("This cell part can not be deleted")
    else:
        for k, v in total_bounds.items():
            for i in range(r-1,r+1):
                for j in range(c-1,c+1):
                    if i == r and j == c:
                        n = str(r) + str(c)
                        if n in k:
                            react[k].bounds = (0.0, 0.0)
                    else:
                        m = str(i) + str(j)
                        if m in k:
                            react[k].bounds = (0.0, v)
    '''
    objective_functions = [v.flux_expression for k, v in react.items() if 'rBcBio' in k]
    quadratic_objective = test_model.problem.Objective(
        sum(objective_functions),
        direction='max')
    test_model.objective = quadratic_objective
    solution = test_model.optimize(objective_sense=None)
    print(solution)
    print({k:v.flux for k,v in react.items()})

    return solution, react
