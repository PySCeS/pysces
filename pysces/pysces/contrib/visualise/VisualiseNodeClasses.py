'''
Colin Gillespie (c.gillespie@ncl.ac.uk)
Last modified: 1/4/05

If want change the shape/colour/etc of a node or arrow
make your changes here
'''
'''
Small changes to format/style and stoichiometry marked as # brett 20050702
Brett G. Olivier (bgoli@users.sourceforge.net)
'''


class EdgeLines:
	'''
	If we want differnt arrow types for the different reaction types
	For example, square 'arrows' for a product
	Just now all the arrows are the same
	'''
	
	def __init__(self, edge_attribute):
		self.edge_attribute = edge_attribute
	
	def reactant(self):
		'''See class doc'''
		if self.edge_attribute['arrowtail'] != 'none':
			self.edge_attribute['arrowtail'] = 'normal'
		self.edge_attribute['arrowhead'] = 'none' # brett 20050702
		self.edge_attribute['arrowsize'] = '0.8'  # brett 20050702
		return self.edge_attribute
	
	def product(self):
		'''See class doc'''
		#if self.edge_attribute['arrowtail'] != 'none':  # brett 20050702
		#	self.edge_attribute['arrowtail'] = 'normal' # brett 20050702
		self.edge_attribute['arrowhead'] = 'normal'

		self.edge_attribute['arrowtail'] = 'none' # brett 20050702
		self.edge_attribute['arrowsize'] = '0.8'  # brett 20050702
		return self.edge_attribute
	
	def modifier(self):
		'''See class doc'''
		#if self.edge_attribute['arrowtail'] != 'none': # brett 20050702
		#	self.edge_attribute['arrowtail'] = 'normal' # brett 20050702
		self.edge_attribute['arrowhead'] = 'normal'
		self.edge_attribute['arrowtail'] = 'none'               # brett 20050702
		self.edge_attribute['weight'] = 0 						# brett 20050702
		self.edge_attribute['style'] = 'dotted, setlinewidth(2)'# brett 20050702
		self.edge_attribute['arrowsize'] = '0.8'				# brett 20050702
		self.edge_attribute['color'] = 'mediumblue' 			# brett 20050702
		return self.edge_attribute

class DotNodes:
	'''A class used to define the node shapes, colours, fonts, etc.
	So you can change what a species node looks like here.
	For edges, see the EdgeLines class'''
	
	def	__init__(self, model):
		self.general_node = {'style':'filled', 'fontsize':10, 'fontname':'Arial'}
		self.m = model
	
	def edgeNode(self, reaction, stoic=0):
		'''See class doc'''
		edge_attribute = {'arrowtail':'none', 'arrowhead':'normal',
			'color':'black', 'fontsize':10, 'fontname':'Arial'}
		
		if reaction.getReversible() == 1:
			edge_attribute['arrowtail'] = 'normal'
		
		if reaction.getFast() == 1: 
			edge_attribute['color'] = 'red'
		
		# if stoic != 1:				# brett 20050702
		if stoic > 0.0 and stoic != 1: # brett 20050702: only show non-int, non-unitary stoich.  
			edge_attribute['label'] = str(stoic) # brett 20050702
		
		return EdgeLines(edge_attribute)

	def specieNode(self, i, other_attributes=None):
		'''See class doc'''
		if not other_attributes:
			other_attributes = {}
		
		label = get_label(self.m.getSpecies(i-1))

		node_attribute = self.general_node
		node_attribute['color'] = 'orange'
		node_attribute['shape'] = 'ellipse'
		node_attribute['label'] = label
		node_attribute['URL'] = "Javascript:species('//sbml:species[%s]')" % str(i)
		
		#if boundary condition
		if self.m.getSpecies(i-1).getBoundaryCondition():
			node_attribute['color'] = "coral"
		
		for key, value in other_attributes.items():
			node_attribute[key] = value
		return node_attribute
	
	def reactionNode(self, i, other_attributes=None):
		'''See class doc'''
		if not other_attributes:
			other_attributes = {}
		
		label = get_label(self.m.getReaction(i-1))
		
		node_attribute = self.general_node
		node_attribute['color'] = 'moccasin'
		node_attribute['shape'] = 'rectangle'
		node_attribute['label'] = label
		node_attribute['URL'] = "Javascript:reaction('//sbml:reaction[%s]')" % str(i)
		
		for key, value in other_attributes.items():
			node_attribute[key] = value
		return node_attribute
		
	def compartmentNode(self, i, other_attributes=None):	
		'''See class doc'''
		if not other_attributes:
			other_attributes = {}
		
		label = get_label(self.m.getCompartment(i-1))
		
		node_attribute = self.general_node
		node_attribute['graph_name'] = self.m.getCompartment(i-1).getId()
		node_attribute['label'] = label
		node_attribute['URL'] = "Javascript:compartment(\'//sbml:compartment[%s]\')" \
			% str(i)
		
		for key, value in other_attributes.items():
			node_attribute[key] = value
		return node_attribute
		
def get_label(obj):
	'''Just now the visualiser displays id's not names
	Uncomment the lines below to display names'''
#	if len(obj.getName()):
#		return obj.getName()
#	else:
#		return obj.getId()
	return obj.getId()
