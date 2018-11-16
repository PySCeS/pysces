'''
Colin Gillespie (c.gillespie@ncl.ac.uk)
Last modified: 4/4/05

The access point to the visualisation stuff
'''

from .VisualiseClasses import VisualiseSpecies, VisualiseReaction, VisualiseModel

class VisualiseFormat:
	'''Static returns a single value
		web returns a piccy with an image map
	'''
	
	formats = ['ps', 'ps2', 'hpgl', 'pcl', 'mif', 'pic', 'gd', 'gd2', 'gif',
		'jpg', 'jpeg', 'png', 'wbmp', 'ismap', 'imap', 'cmap', 'vrml', 'vtx', 
		'mp','fig', 'svg', 'svgz', 'dia', 'dot', 'canon', 'plain', 'plain-ext',
		'xdot']
	
	def __init__(self, py_dot_graph):
		self.pyGraph = py_dot_graph
	
	def web(self, frmt='gif'):
		'''Returns a base64 encoded graph and an image map'''
		if frmt not in self.formats:
			raise ValueError('Format not recognised.')
		graph = self.pyGraph.create(format=frmt)
		cmap = self.pyGraph.create(format='cmap')
		return graph, cmap
	
	def static(self, frmt='dot'):
		'''returns just one value'''
		if frmt not in self.formats:
			raise ValueError('Format not recognised.')
		return self.pyGraph.create(format=frmt)
		
class VisualiseSBML:
	"""The main visualisation class"""
	
	def	__init__(self, m):
		self.m = m
		
	def visualiseModel(self):
		'''Visualise the entire SBML model'''
		py_dot_graph = VisualiseModel(self.m).getPyDot()
		return VisualiseFormat(py_dot_graph)
		
	def visualiseReaction(self, rea_num):
		'''
		Visualise a single reaction
		rea_num is the number of the reaction you want to visualise
		Numbering starts at 1
		'''
		py_dot_graph = VisualiseReaction(self.m, rea_num).getPyDot()
		return VisualiseFormat(py_dot_graph)

	def visualiseSpecies(self, sp_num):
		'''
		Visualise a single species
		sp_num is the number of the species you want to visualise
		Numbering starts at 1
		'''
		py_dot_graph = VisualiseSpecies(self.m, sp_num).getPyDot()
		return VisualiseFormat(py_dot_graph)

if __name__ == "__main__":
	sbml_file = open('../chapmodel1.xml','r').read()
	from BasisTools import BasisTools
	sbml1 = BasisTools(sbml = sbml_file)
	sbml =  sbml1.m

	print('model visualise')
	print(VisualiseSBML(sbml).visualiseModel().static('dot'))

	print('Specie visualise')
	print(VisualiseSBML(sbml).visualiseSpecies(1).static('dot'))

	print('Reaction visualise')
	print(VisualiseSBML(sbml).visualiseReaction(2).static('dot'))
	
	print('Empty Model')
	sbml1 = BasisTools()
	sbml = sbml1.m
	print(VisualiseSBML(sbml).visualiseModel().static('dot'))
	
	print('One Compartment 1 Species')
	sbml1 = BasisTools()
	sbml1.add(Compartments={'Id':'Cell'})
	sbml1.add(Species={'Id':'sp1', 'Compartment':'Cell'})
	sbml = sbml1.m
	print(VisualiseSBML(sbml).visualiseModel().static('dot'))
	
	sbml =  sbml1.m
	print('EOF')
	
