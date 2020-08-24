"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2020 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""

from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from .version import __version__

__doc__ = '''PySCeS ModelMap module: useful for exploring model component relations'''


class ModelMapBase(object):
    name = None

    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name

    def get(self, attr):
        """Return an attribute whose name is str(attr)"""
        try:
            return getattr(self, attr)
        except:
            print("%s is not an attribute of this instance" % attr)
            return None


class MapList(list):
    def __init__(self, *args):
        list.__init__(self, *args)

    def asSet(self):
        return set(self.__getslice__(0, self.__len__()))


class ModelMap(ModelMapBase):
    __nDict__ = None
    reactions = None
    species = None
    species_fixed = None
    compartments = None
    __model__ = None
    __InitStrings__ = None
    __InitDict__ = None
    __not_inited__ = None
    global_parameters = None
    __parameter_store__ = None

    def __init__(self, model):
        self.setName(model.ModelFile[:-4])
        self.__nDict__ = model.__nDict__
        self.__model__ = model
        self.__InitDict__ = self.__model__.__InitDict__.copy()
        self.__compartments__ = self.__model__.__compartments__.copy()
        for k in list(self.__InitDict__.keys()):
            self.__InitDict__[k] = getattr(self.__model__, k)
        self.global_parameters = []
        self.__parameter_store__ = []
        self.__not_inited__ = []

        # operational shortcuts
        self.addSpecies()
        self.addReactions()
        self.generateMappings()
        self.addCompartments()

    def __cleanString__(self, s):
        s = s.lstrip()
        s = s.rstrip()
        return s

    def addOneSpecies(self, species, fix=False):
        s = Species(species)
        s.setValue(self.__InitDict__[s.name])
        if fix:
            s.fixed = True
        setattr(self, species, s)
        self.species.append(s)
        if fix:
            self.species_fixed.append(s)

    def addCompartments(self):
        self.compartments = []
        for c in self.__model__.__compartments__:
            co = Compartment(
                self.__model__.__compartments__[c]['name'],
                self.__model__.__compartments__[c]['size'],
            )
            self.compartments.append(co)
            setattr(self, c, co)
        cname = [c.name for c in self.compartments]
        for s in list(self.__model__.__sDict__.keys()):
            if self.__model__.__sDict__[s]['compartment'] in cname:
                getattr(self, self.__model__.__sDict__[s]['compartment']).setComponent(
                    getattr(self, s)
                )
                getattr(self, s).compartment = getattr(
                    self, self.__model__.__sDict__[s]['compartment']
                )
        for r in list(self.__model__.__nDict__.keys()):
            if self.__model__.__nDict__[r]['compartment'] in cname:
                getattr(self, self.__model__.__nDict__[r]['compartment']).setComponent(
                    getattr(self, r)
                )
                getattr(self, r).compartment = getattr(
                    self, self.__model__.__nDict__[r]['compartment']
                )

    def addOneReaction(self, reaction):
        r = Reaction(reaction)
        r.addFormula(self.__nDict__[r.name]['RateEq'].replace('self.', ''))
        if self.__nDict__[r.name]['Type'] == 'Irrev':
            r.reversible = False

        fxnames = self.hasFixedSpecies()
        for p in self.__nDict__[r.name]['Params']:
            p = p.replace('self.', '')
            if (
                p not in self.hasGlobalParameters()
                and p not in fxnames
                and p not in self.__compartments__
            ):
                if p in self.__InitDict__:
                    par = Parameter(p, self.__InitDict__[p])
                else:
                    par = Parameter(p)
                    if p not in self.__not_inited__:
                        self.__not_inited__.append(p)
                par.setAssociation(r)
                self.global_parameters.append(par)
                setattr(self, p, par)
                r.addParameter(par)
            elif p not in fxnames and p not in self.__compartments__:
                pidx = self.hasGlobalParameters().index(p)
                self.global_parameters[pidx].setAssociation(r)
                r.addParameter(self.global_parameters[pidx])
        setattr(self, reaction, r)
        self.reactions.append(r)

    def addSpecies(self):
        self.species = []
        self.species_fixed = []
        for s in self.__model__.species:
            self.addOneSpecies(s, fix=False)
        for s in self.__model__.fixed_species:
            self.addOneSpecies(s, fix=True)

    def addReactions(self):
        self.reactions = []
        for r in self.__model__.reactions:
            self.addOneReaction(r)

    def generateMappings(self):
        for reac in self.reactions:
            for reag in self.__nDict__[reac.name]['Reagents']:
                if self.__nDict__[reac.name]['Reagents'][reag] < 0.0:
                    reac.addSubstrate(getattr(self, reag.replace('self.', '')))
                    getattr(self, reag.replace('self.', '')).setSubstrate(
                        getattr(self, reac.name)
                    )
                else:
                    reac.addProduct(getattr(self, reag.replace('self.', '')))
                    getattr(self, reag.replace('self.', '')).setProduct(
                        getattr(self, reac.name)
                    )
                reac.stoichiometry.setdefault(
                    reag.replace('self.', ''),
                    self.__nDict__[reac.name]['Reagents'][reag],
                )
            for mod in self.__nDict__[reac.name]['Modifiers']:
                reac.addModifier(getattr(self, mod.replace('self.', '')))
                getattr(self, mod.replace('self.', '')).setModifier(
                    getattr(self, reac.name)
                )

    def hasReactions(self):
        return MapList([r.name for r in self.reactions])

    def hasSpecies(self):
        return MapList([s.name for s in self.species])

    def hasFixedSpecies(self):
        return MapList([s.name for s in self.species_fixed])

    def findReactionsThatIncludeAllSpecifiedReagents(self, *args):
        assert len(args) > 1, '\nNeed two or more species for this one!'
        setlist = [getattr(self, s).isReagentOf().asSet() for s in args]
        isect = setlist[0]
        for s in setlist:
            isect.intersection_update(s)
        return MapList(isect)

    def hasGlobalParameters(self):
        return MapList(p.name for p in self.global_parameters)


class Reaction(ModelMapBase):
    modifiers = None
    substrates = None
    products = None
    stoichiometry = None
    parameters = None
    reversible = True
    formula = None
    compartment = None

    def __init__(self, name):
        self.setName(name)
        self.modifiers = []
        self.substrates = []
        self.products = []
        self.stoichiometry = {}
        self.parameters = []

    def addSubstrate(self, species):
        setattr(self, species.name, species)
        self.substrates.append(species)

    def addProduct(self, species):
        setattr(self, species.name, species)
        self.products.append(species)

    def addModifier(self, species):
        setattr(self, species.name, species)
        self.modifiers.append(species)

    def addFormula(self, formula):
        self.formula = formula

    def addParameter(self, par):
        setattr(self, par.name, par)
        self.parameters.append(par)

    def hasProducts(self, t=type):
        return MapList([p.name for p in self.products])

    def hasSubstrates(self):
        return MapList([s.name for s in self.substrates])

    def hasModifiers(self):
        return MapList([m.name for m in self.modifiers])

    def hasParameters(self):
        return MapList([p.name for p in self.parameters])

    def hasReagents(self):
        return MapList(self.hasSubstrates() + self.hasProducts())


class NumberBase(ModelMapBase):
    value = None

    def __call__(self):
        return self.value

    def getValue(self):
        return self.value

    def setValue(self, v):
        self.value = v


class Species(NumberBase):
    subs = None
    prods = None
    mods = None
    fixed = False
    compartment = None

    def __init__(self, name):
        self.setName(name)
        self.subs = []
        self.prods = []
        self.mods = []

    def setSubstrate(self, reaction):
        setattr(self, reaction.name, reaction)
        self.subs.append(reaction)

    def setProduct(self, reaction):
        setattr(self, reaction.name, reaction)
        self.prods.append(reaction)

    def setModifier(self, reaction):
        setattr(self, reaction.name, reaction)
        self.mods.append(reaction)

    def isSubstrateOf(self):
        return MapList([r.name for r in self.subs])

    def isProductOf(self):
        return MapList([r.name for r in self.prods])

    def isModifierOf(self):
        return MapList([r.name for r in self.mods])

    def isReagentOf(self):
        return MapList(self.isSubstrateOf() + self.isProductOf())


class Parameter(NumberBase):
    association = None
    formula = None

    def __init__(self, name, value=None):
        self.name = name
        self.value = value
        self.association = []

    def setAssociation(self, reac):
        self.association.append(reac)
        setattr(self, reac.name, reac)

    def isParameterOf(self):
        return MapList([a.name for a in self.association])

    def setFormula(self, formula):
        self.formula = formula


class Compartment(NumberBase):
    components = None

    def __init__(self, name, value=None):
        self.name = name
        self.value = value
        self.components = []

    def setComponent(self, comp):
        self.components.append(comp)
        setattr(self, comp.name, comp)

    def hasComponents(self):
        return MapList([a.name for a in self.components])


if __name__ == '__main__':
    import pysces

    M = pysces.model('pysces_model_linear1')
    M.doLoad()

    print('\nModel', M.ModelFile)
    print('=============')
    modmap = ModelMap(M)

    print('Reactions\n', modmap.hasReactions())
    print('Species\n', modmap.hasSpecies())
    print('FixedSpecies\n', modmap.hasFixedSpecies())
    print(' ')
    print('R1 has reagents\n', modmap.R1.hasReagents())
    print('R1 has sub\n', modmap.R1.hasSubstrates())
    print('R1 has prod\n', modmap.R1.hasProducts())
    print('R1 has mod\n', modmap.R1.hasModifiers())
    print(' ')
    print('s2 is reagent\n', modmap.s2.isReagentOf())
    print('s2 is sub\n', modmap.s2.isSubstrateOf())
    print('s2 is prod\n', modmap.s2.isProductOf())
    print('s2 is mod\n', modmap.s2.isModifierOf())
    print(' ')
    print('R2 stoich\n', modmap.R2.stoichiometry)
    print(' ')
    print(
        'findReactionsThatIncludeAllSpecifiedReagents(A, B):',
        modmap.findReactionsThatIncludeAllSpecifiedReagents('s1', 's2'),
    )

    print('\nmodmap.hasGlobalParameters\n', modmap.hasGlobalParameters())

    print('\nParameter associations')
    for p in modmap.global_parameters:
        print('%s.isParameterOf() %s' % (p.name, p.isParameterOf()))
