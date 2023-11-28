""" Methods for running discretised FBA

:Author: Yin Hoon Chew <yinhoon.chew@bham.ac.uk>
:Author: Alexandra Anamaria-Sorinca
:Date: 13-09-2022
:License: MIT

"""

from itertools import product, starmap
import cobra
import numpy as np
import random


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
            aspect_ratio (:obj:`float`): length divided by width
            peri_to_area (:obj:`float`): number of regions in the outer
                layer divided by the total number of regions
            regions (:obj:`numpy.object`): array of numpy objects with a
                size of length x width
            model (:obj:`cobra.Model`): FBA model object
            cell_level_reactions (:obj:`dict`): dictionary of cell-level 
                reactions, i.e. reactions in the extracellular compartment         

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
        self.aspect_ratio = self.length / self.width
        
        if self.width==1 or self.length==1:
            pa_ratio = 1
        else:    
            pa_ratio = (2*(self.width + self.length) - 4) / (self.width*self.length)
        self.peri_to_area = pa_ratio
        
        self.regions = np.empty((self.length, self.width), dtype=object)
        self.model = cobra.Model()
        self.cell_level_reactions = {}
        for i in range(self.length):
            for j in range(self.width):
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
        for m in extra_mets:
            metabolite_objs.append(cobra.Metabolite(m))
        self.model.add_metabolites(metabolite_objs)           
        
        # create reactions in the extracellular compartment, i.e. cell-level reactions
        reaction_objs = []
        for k in extra_rxns.keys():
            rxn_obj = cobra.Reaction(k)
            reaction_objs.append(rxn_obj)
        self.model.add_reactions(reaction_objs)    
        for rxn_obj in reaction_objs:
            rxn_obj.reaction = extra_rxns[rxn_obj.id]
            self.cell_level_reactions[rxn_obj.id] = rxn_obj                    
        
        # create intracellular metabolites for each discrete region
        metabolite_objs = []
        for i in range(self.length):
            for j in range(self.width):
                region = self.regions[i, j]
                for m in intra_mets:
                    region_metabolite = cobra.Metabolite(f'{m}_{region.id}')
                    metabolite_objs.append(region_metabolite)
                    region.add_metabolite(m, region_metabolite)    
        self.model.add_metabolites(metabolite_objs)

        # create intracellular reactions for each discrete region
        for i in range(self.length):
            for j in range(self.width):
                region = self.regions[i, j]
                for k, v in intra_rxns.items():
                    region_rxn = cobra.Reaction(f'{k}_{region.id}')                    
                    region.add_synthesis_reaction(k, region_rxn)
                    self.model.add_reactions([region_rxn])
                    new_rxn_str = v 
                    for m in intra_mets:
                        if m in new_rxn_str:
                            new_rxn_str = new_rxn_str.replace(m, f'{m}_{region.id}')
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
        x,y = np.ogrid[0:self.length, 0:self.width]
        # retrieve the inner layers
        inner_layers = (x>0)&(x<self.length-1)&(y>0)&(y<self.width-1)
        # create transport reaction into the outermost layer
        for region in self.regions[~inner_layers]:
            for k, v in extra_intra_rxns.items():
                region_rxn = cobra.Reaction(f'{k}_{region.id}')                    
                region.add_transport_reaction(k, region_rxn)
                self.model.add_reactions([region_rxn])
                new_rxn_str = v 
                for m in intra_mets:
                    if m in new_rxn_str:
                        new_rxn_str = new_rxn_str.replace(m, f'{m}_{region.id}')
                region_rxn.reaction = new_rxn_str
        
    def create_diffusion(self, intra_mets):
        """ Add nutrient diffusion between neighbouring regions

        Args:
            intra_mets (:obj:`list` of :obj:`str`): list of IDs of intracellular 
                metabolites that diffuse in the cell
        """
        for i in range(self.length):
            for j in range(self.width):
                region = self.regions[i, j]
                neighbours = self.get_neighbouring_regions(
                    region.row_number, region.column_number)
                for neighbour in neighbours:
                    neighbour_region = self.regions[neighbour[0], neighbour[1]]
                    for m in intra_mets:            
                        region_rxn = cobra.Reaction(
                            f'{m}_{region.id}_to_{neighbour_region.id}')                    
                        region.add_diffusion(region_rxn)
                        neighbour_region.add_diffusion(region_rxn)
                        self.model.add_reactions([region_rxn])
                        region_rxn.reaction = \
                            f'{m}_{region.id} --> {m}_{neighbour_region.id}'                

    def get_neighbouring_regions(self, row_index, col_index):
        """ Find all the neighbours of a region

        Args:
            row_index (:obj:`int`): row index of the region
            col_index (:obj:`int`): column index of the region

        Returns:
            :obj:`list`: list of row and column indices of all neighbours     
        """         
        length, width = len(self.regions), len(self.regions[0])
        neighbours = list(starmap(lambda a, b: (row_index + a, col_index + b), product((0, -1, +1), (0, -1, +1))))
        neighbours.pop(0) #  exclude region of interest
        neighbours = list(filter(lambda cell: cell[0] in range(length) and cell[1] in range(width), neighbours))
        return neighbours

    def distribute_enzyme(self, enzyme_id, total_concentration, gradient=0, 
        random_distribution=False, random_seed=None, all_regions=True):
        """ Distribute an enzyme into each region based on the distribution type

        Args:
            enzyme_id (:obj:`str`): enzyme ID
            total_concentration (:obj:`float`): cellular concentration of the enzyme
            gradient (:obj:`float`, optional): the gradient of enzyme distribution;
                only takes a value between -1.0 to 1.0; a higher absolute value will 
                create a steeper gradient; a positive value will set the concentration 
                to increase from outer to inner layers; a negative value will set the 
                concentration to decrease from outer to inner layers; a zero value will 
                distribute the concentration uniformly (default)
            random_distribution (:obj:`bool`, optional): if True, enzyme will be randomly
                distributed; default is False
            random_seed (:obj:`int`, optional): the seed value used for random number 
                generator if random distribution is true. If no value is provided,
                tha random number will not be seeded    
            all_regions (:obj:`bool`, optional): if True, enzyme is distributed to all
                regions, else enzyme is only distributed to regions in the outermost 
                layer; default is True
        
        Raises:
            :obj:`ValueError`: if the value of gradient is outside the range between 
                -1.0 and 1.0
        """
        if gradient < -1 or gradient > 1:
            raise ValueError('gradient can only take a value between -1.0 and 1.0')
        
        if all_regions:

            if random_distribution:
                if random_seed:
                    random.seed(random_seed)
                sample = [random.randint(0,100) for i in range(self.width*self.length)]
                enzyme_distribution = [total_concentration*i/sum(sample) for i in sample]
                count = 0
                for i in range(self.length):
                    for j in range(self.width):
                        self.regions[i, j].enzymes[enzyme_id] = enzyme_distribution[count]
                        count += 1
            
            else:
                if gradient == 0:
                    distribution = np.ones((self.length, self.width), dtype=float)                    
                
                else:
                    if gradient > 0:
                        start_value = 1
                    else:    
                        start_value = min(self.length, self.width) / 2

                    distribution = np.full((self.length, self.width), start_value, 
                        dtype=float)
                    x,y = np.ogrid[0:self.length, 0:self.width]
                    z = 1
                    while z < self.length and z < self.width:
                        inner_layers = (x>0+z-1)&(x<self.length-z)&(y>0+z-1)&(
                            y<self.width-z)
                        for i in range(self.length):
                            for j in range(self.width):
                                if inner_layers[i, j]:
                                    distribution[i, j] += gradient                    
                        z += 1
                for i in range(self.length):
                    for j in range(self.width):     
                        self.regions[i, j].enzymes[enzyme_id] = \
                            total_concentration*distribution[i, j]/sum(sum(distribution))
        
        else:
            x,y = np.ogrid[0:self.length, 0:self.width]
            inner_layers = (x>0)&(x<self.length-1)&(y>0)&(y<self.width-1)
            outermost_layer = self.regions[~inner_layers]            
            if random_distribution:
                if random_seed:
                    random.seed(random_seed)
                sample = [random.randint(0,100) for i in range(len(outermost_layer))]
                enzyme_distribution = [total_concentration*i/sum(sample) for i in sample]
                count = 0
                for region in outermost_layer:
                    region.enzymes[enzyme_id] = enzyme_distribution[count]
                    count += 1
            else:
                uniform_concentration = total_concentration/len(outermost_layer)
                for region in outermost_layer:
                    region.enzymes[enzyme_id] = uniform_concentration
            for region in self.regions[inner_layers]:
                region.enzymes[enzyme_id] = 0            

    def calc_enzymatic_bounds(self, enzyme_rxn_assoc):
        """ Calculate the upper bound of each reaction in each region
            using the formula k_cat * [E] where k_cat is the catalytic
            constant and [E] is the concentration of enzyme. If there are
            more than one enzyme catalysing a reaction, the bound will
            the sum of all.

        Args:
            enzyme_rxn_assoc (:obj:`dict` of :obj:`dict`): the dictionary of
                enzyme (ID) and its rate constant that catalyse a reaction (ID as key)

        Returns:
            :obj:`dict`: the upper bound of each enzyme catalysed synthesis reaction
            :obj:`dict`: the upper bound of each enzyme catalysed transport reaction        
        """
        synthesis_reaction_bounds = {}
        transport_reaction_bounds = {}
        for i in range(self.length):
            for j in range(self.width):
                region = self.regions[i, j]
                for k, v in region.synthesis_reactions.items():
                    if k in enzyme_rxn_assoc:
                        upper_bound = sum([k_cat*region.enzymes[enz] \
                            for enz, k_cat in enzyme_rxn_assoc[k].items()])
                        synthesis_reaction_bounds[v.id] = upper_bound
                for k, v in region.transport_reactions.items():
                    if k in enzyme_rxn_assoc:
                        upper_bound = sum([k_cat*region.enzymes[enz] \
                            for enz, k_cat in enzyme_rxn_assoc[k].items()])
                        transport_reaction_bounds[v.id] = upper_bound

        return synthesis_reaction_bounds, transport_reaction_bounds                                

    def set_bounds(self, bounds):
        """ Set the bounds for all reactions and diffusions

        Args:
            bounds (:obj:`dict`): tuple of max and min bounds for each given
                reaction/diffusion ID
        """
        for rxn, bounds in bounds.items():
            self.model.reactions.get_by_id(rxn).bounds = bounds
                
    def solve(self):
        """ Optimise the objective function

        Returns:
            :obj:`cobra.core.solution`: Optimised solution
        """
        return self.model.optimize(objective_sense=None)    
            

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
        self.synthesis_reactions = {}
        self.transport_reactions = {}
        self.diffusions = {}
        self.enzymes = {}

    def add_metabolite(self, metabolite_id, metabolite):
        """ Add a metabolite to the region

        Args:
            metabolite_id (:obj:`str`): metabolite ID
            metabolite (:obj:`cobra.Metabolite`): a metabolite object
        """
        self.metabolites[metabolite_id] = metabolite

    def add_synthesis_reaction(self, reaction_id, reaction):
        """ Add a synthesis reaction to the region

        Args:
            reaction_id (:obj:`str`): reaction ID
            reaction (:obj:`cobra.Reaction`): a reaction object
        """
        self.synthesis_reactions[reaction_id] = reaction

    def add_transport_reaction(self, reaction_id, reaction):
        """ Add a synthesis reaction to the region

        Args:
            reaction_id (:obj:`str`): reaction ID
            reaction (:obj:`cobra.Reaction`): a reaction object
        """
        self.transport_reactions[reaction_id] = reaction    

    def add_diffusion(self, diffusion):
        """ Add diffusion of a metabolite to the region

        Args:
            diffusion (:obj:`cobra.Reaction`): a reaction object
        """
        self.diffusions[diffusion.id] = diffusion

    def add_enzyme_concentration(self, enzyme_id, concentration):
        """ Add an enzyme and its concentration to the region

        Args:
            enzyme_id (:obj:`str`): enzyme ID
            concentration (:obj:`float`): concentration of enzyme in the region
        """
        self.enzymes[enzyme_id] = concentration 
