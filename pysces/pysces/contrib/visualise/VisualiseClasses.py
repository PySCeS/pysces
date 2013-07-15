'''
Colin Gillespie (c.gillespie@ncl.ac.uk)
Last modified: 1/4/05

If want change the shape/colour/etc of a node or arrow
look at VisualiseNodeClasses
'''
'''
Small changes to stoichiometry and modifiers marked as # brett 20050702
Brett G. Olivier (bgoli@users.sourceforge.net)
'''
import pydot
from VisualiseNodeClasses import DotNodes, EdgeLines

class VisualiseModel(DotNodes, EdgeLines):
	"""The visualisation model class"""
	
	def	__init__(self, m):
		'''
		m is a libsbml getModel file
		initialise a few short-cuts here
		'''
		self.m = m
		self.py_dot_graph = pydot.Dot(type='digraph', graph_name='graphName', \
			rankdir = "LR")
		self.add_edge = self.py_dot_graph.add_edge
		self.add_node = self.py_dot_graph.add_node
		DotNodes.__init__(self, m)

	def __createCompartmentDict(self):
		"""creates a dictionary with all compartments as keys and a 
			pydot cluster as a value
		"""
		compart_dict = {}
		for i, co in enumerate(self.m.getListOfCompartments()):
			compart_dict[co.getId()] = [co.getOutside()]
			node_attribute = self.compartmentNode(i+1)
			compart_dict[co.getId()].append(pydot.Cluster(**node_attribute))
		return compart_dict
	
	def __organiseCompartments(self, compart_dict):
		"""organises the compartments, using outside"""
		for key in compart_dict.keys():
			outside = compart_dict[key][0]
			if len(outside) ==0:
				self.py_dot_graph.add_subgraph(compart_dict[key][1])
			else:
				compart_dict[outside][1].add_subgraph(compart_dict[key][1])
		return self.py_dot_graph	

	def __isGraph2Big(self):
		'''Is the model too big'''
		if (self.m.getNumReactions() + self.m.getNumSpecies()) < 250:
			return 0

		dict1 = {'URL':'#', 'label':'There are too many reactions/ \
					species to be display'}
		node_attribute = self.compartmentNode(1, dict1)
		cluster = pydot.Cluster(**node_attribute)
		
		node_attribute = self.reactionNode(2, {'URL':'#'})
		node_attribute['label'] = '%d Reactions' % self.m.getNumReactions()
		cluster.add_node(pydot.Node('reactions', **node_attribute))
		
		node_attribute = self.specieNode(2, {'URL':'#'})
		node_attribute['label'] = '%d species' % self.m.getNumSpecies()
		cluster.add_node(pydot.Node('species', **node_attribute))
		
		self.py_dot_graph.add_subgraph(cluster)
		return 1
			
	def __isCompartments(self):
		'''Does the model have any compartments'''
		if self.m.getNumCompartments() > 0:
			return 0
		
		node_attribute = {'graph_name':'empty model', 'label':'empty model', \
				'fontsize':10, 'fontname':'Arial'}
		node_attribute['URL'] = '#'
		cluster = pydot.Cluster(**node_attribute)
			
		node_attribute = {'color':'moccasin', 'shape':'rectangle',
			'style':'filled', 'label':'0 species', 'fontsize':10}
		node_attribute['fontname'] = 'Arial'
		node_attribute['URL'] = '#'
		cluster.add_node(pydot.Node('species', **node_attribute))
	
		self.py_dot_graph.add_subgraph(cluster)
		return 1
	
	def __isSpecies(self):
		'''Does the model have any species'''
		if self.m.getNumSpecies() > 0:
			return 0
		
		node_attribute = self.compartmentNode(1)
		cluster = pydot.Cluster(**node_attribute)
		node_attribute = {'color':'moccasin', 'shape':'rectangle',
			'style':'filled', 'label':'0 species', 'fontsize':10}
		node_attribute['fontname'] = 'Arial'
		node_attribute['URL'] = '#'
		cluster.add_node(pydot.Node('species', **node_attribute))
		
		self.py_dot_graph.add_subgraph(cluster)
		return 1
	
	def __addReactionToCompartment(self, sp_ref):
		'''Have to put the reaction into a compartment so
		choose the first reactant and stick it that compartment
		If there are no reactants than get the first product
		'''
		for sp in self.m.getListOfSpecies():
			if sp.getId() == sp_ref.getSpecies():
				return 1
	
	def getPyDot(self):
		'''main method'''
		#If there are too many reactions
		
		if self.__isGraph2Big():
			return self.py_dot_graph
		
		if self.__isCompartments():
			return self.py_dot_graph
			
		if self.__isSpecies():
			return self.py_dot_graph

		

		compart_dict = self.__createCompartmentDict()
		#Add the species into a compartment
		for i, sp in enumerate(self.m.getListOfSpecies()):
			node_attribute = self.specieNode(i+1)
			compart_dict[sp.getCompartment()][1].add_node(
				pydot.Node(sp.getId(), **node_attribute))
	
		for i, reaction in enumerate(self.m.getListOfReactions()):
			node_attribute = self.reactionNode(i+1)
			for j, sp_ref in enumerate(reaction.getListOfReactants()):
				if j == 0 and self.__addReactionToCompartment(sp_ref):
					compart_dict[sp.getCompartment()][1].add_node(
						pydot.Node(reaction.getId(), **node_attribute))
				
				# stoc = int(sp_ref.getStoichiometry()) # brett 20050702
				stoc = sp_ref.getStoichiometry() 		# brett 20050702: allow for non-int stoich
				
				edge_attribute = self.edgeNode(reaction, stoc).reactant()
				edge = pydot.Edge(sp_ref.getSpecies(), 
					reaction.getId(), **edge_attribute)
				self.add_edge(edge)
				
			for j, sp_ref in enumerate(reaction.getListOfProducts()):
				if j == 0 and not reaction.getNumReactants() and \
					self.__addReactionToCompartment(sp_ref):
					compart_dict[sp.getCompartment()][1].add_node(
							pydot.Node(reaction.getId(), **node_attribute))
						
				#stoc = int(sp_ref.getStoichiometry()) 	# brett 20050702
				stoc = sp_ref.getStoichiometry() 		# brett 20050702: allow for non-int stoich
			
				edge_attribute = self.edgeNode(reaction, stoc).product()
				edge = pydot.Edge(reaction.getId(), 
						sp_ref.getSpecies(), **edge_attribute)
				self.add_edge(edge)
	
			for sp_ref in reaction.getListOfModifiers():
				edge_attribute = self.edgeNode(reaction).modifier()
				edge = pydot.Edge(sp_ref.getSpecies(), reaction.getId(), **edge_attribute)
				self.add_edge(edge)
				# perhaps modifiers should only go from species to reaction # brett 20050702
				# edge = pydot.Edge(reaction.getId(), sp_ref.getSpecies(), **edge_attribute) # brett 20050702
				# self.add_edge(edge) # brett 20050702

		return self.__organiseCompartments(compart_dict)
	
class VisualiseReaction(DotNodes, EdgeLines):
	'''Used to visualise a single reaction'''
	
	def	__init__(self, m, re_num):
		'''
		m is libsbml.getModel 
		re_num is the reaction to be visualised
		'''
		
		self.m = m
		self.py_dot_graph = pydot.Dot(type='digraph', graph_name='graphName', \
			rankdir = "LR")
		self.add_edge = self.py_dot_graph.add_edge
		self.add_node = self.py_dot_graph.add_node
		DotNodes.__init__(self, m)
		
		self.re_num = re_num
		self.re = self.m.getReaction(re_num-1)
		
		DotNodes.__init__(self, m)

	def getPyDot(self):
		'''Returns a pyDot object'''
				
		#Create the reaction node
		node_attribute = self.reactionNode(self.re_num)
		self.add_node(pydot.Node(self.re.getId(), **node_attribute))
			
		types = ['Reactants', 'Products', 'Modifiers']
		for typ in types:
			for sp_ref in getattr(self.re, 'getListOf' + typ)():
				sp_ref_id = sp_ref.getSpecies()
				self.__findMatchingSpecies(sp_ref_id)
			
				#stoic = int(sp_ref.getStoichiometry())	# brett 20050702
				stoic = sp_ref.getStoichiometry()  		# brett 20050702: allow for non-int stoich
				
				#XXX: Should be .type not .reactant
				edge_attribute = self.edgeNode(self.re, stoic).reactant()
				if typ != 'Products':
					edge = pydot.Edge(sp_ref_id, self.re.getId(), **edge_attribute)
					self.add_edge(edge)
				if typ != 'Reactants':
					edge = pydot.Edge(self.re.getId(), sp_ref_id, **edge_attribute)
					self.add_edge(edge)
		return self.py_dot_graph			
				
	def __findMatchingSpecies(self, sp_ref_id):
		'''determines where the sp_ref is in listOfSpecies'''
		
		#TODO: make a species list, so we don't keep cycling through sp!
		for i, sp in enumerate(self.m.getListOfSpecies()):
			if sp.getId() == sp_ref_id:
				node_attribute = self.specieNode(i+1)
				self.add_node(pydot.Node(sp_ref_id, **node_attribute))
				break
		return 1

class VisualiseSpecies(DotNodes, EdgeLines):
	'''Visualise a single species'''
	
	def	__init__(self, m, sp_num):
		'''takes in a sbml and validates it'''
		self.m = m
		self.py_dot_graph = pydot.Dot(type='digraph', graph_name='graphName', \
			rankdir = "LR")
		self.add_edge = self.py_dot_graph.add_edge
		self.add_node = self.py_dot_graph.add_node
		DotNodes.__init__(self, m)
		
		self.sp_num = sp_num
		self.sp_id = self.m.getSpecies(sp_num-1).getId()
	
	def getPyDot(self):
		'''If you want to visualise a species'''
		#Create the species node
		node_attribute = self.specieNode(self.sp_num)
		self.add_node(pydot.Node(self.sp_id, **node_attribute))
				
		for i, reaction in enumerate(self.m.getListOfReactions()):
			self.__checkReactionTypesForSpecies(reaction, i)
		return self.py_dot_graph
			
	def __checkReactionTypesForSpecies(self, reaction, reaction_num):
		'''Checks all the species references to determine if the sp is 
		in this reaction
		'''
		types = ['Reactants', 'Products', 'Modifiers']
		for typ in types:
			for sp_ref in getattr(reaction, 'getListOf' + typ)():
				if sp_ref.getSpecies() == self.sp_id:
					node_attribute = self.reactionNode(reaction_num+1)
					r_id = reaction.getId()
					self.add_node(pydot.Node(r_id, **node_attribute))
					
					# stoic = int(sp_ref.getStoichiometry())# brett 20050702
					stoic = sp_ref.getStoichiometry()  		# brett 20050702: allow for non-int stoich
					
					#XXX: Change .modifier
					edge_attribute = self.edgeNode(reaction, stoic).reactant()
					if typ != 'Reactants':
						self.add_edge(pydot.Edge(r_id, self.sp_id, **edge_attribute))
					if typ != 'Products':
						self.add_edge(pydot.Edge(self.sp_id, r_id, **edge_attribute))
					return 1


	
