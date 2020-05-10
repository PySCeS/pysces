import pysces
try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass

# Creates a ratechar object for the lin5_hill model
rc = pysces.RateChar('lin5_hill')

# This function does parameter scans for each species in the patway
rc.doAllRateChar()

# getRateCharData returns the ratechar data for a certain species (C in
# this case). It has some legacy functions and a function to write the 
# data to file and a function to find determine if the species is 
# regulatory for any reactions
rcdC = rc.getRateCharData('C')


print('\n\n\n')
print('C is regulatory for the following reactions:')
# This loop prints each reaction for which C is a regulatory metabolite
for reaction in rcdC.findRegulatory():
    print(reaction)
    print('\n')
    


input('Press Enter to continue...')


# The figure function returns a figure object that is used to plot the 
# data stored in the ratechardata object
figC = rcdC.figure()


# The functions below create specific plots from the data. The 'indiv' 
# functions are useful when there are multiple supply or demand reactions
# as they plot each reaction on it's own figure. 
figC.plotRateChar()
input('Press Enter to continue...')
figC.plotElasRC()
input('Press Enter to continue...')
figC.plotElasRCIndiv()
input('Press Enter to continue...')
figC.plotPartialRCDemand()
input('Press Enter to continue...')
figC.plotPartialRCSupply()
input('Press Enter to continue...')

# The figure object is even more powerful as it also allows you to create
# custom plots. 

# individual lines can  be switched on or off using the set_visible method
# use the reactions names as set out in the reactions dictionary (fig.reactions 
# in this example)
# e for elasticity
# f for flux
# r for response coefficient
# p for partial response coefficients
figC._clear()
# clears the figure from any previous plots

figC.title_label = 'Model: lin5_hill - Green line intersecting a red line'
# changing the value of this variable sets the title of the figure

figC.set_visible('J_R3',f=True)
figC.set_visible('J_R4',f=True)

figC.show()
#shows the figure

figC.save('testfigure1.png')  
#saves the figure with the filename provided in the pysces directory

input('Press ENTER to continue...')



# set_type_visible sets the visibility of all lines of a certain type
# The options are:
#    supply_flux
#    demand_flux
#    supply_rc
#    demand_rc
#    supply_elas
#    demand_elas
#    mod_elas
#    supply_partial_rc
#    demand_partial_rc
figC._clear()
figC.title_label = ''
figC.set_type_visible('supply_flux',True)
figC.set_type_visible('demand_flux',True)
figC.set_type_visible('supply_partial_rc',True,affzero=True) 
## when affzero is true, partial response coefficients with slopes with
## a gradient of zero will also be affected 
figC.show()
figC.save('testfigure2.png') 
input('Press ENTER to continue...')




# In addition to showing the default plots on the screen, we can also save
# them in the ratechar_out directory under the name of the model and with
# subdirectories for each species.
for species in rc._psc_mod.species:
    rcdX = rc.getRateCharData(species)
    figX = rcdX.figure()
    
    figX.plotRateChar(1)
    figX.plotElasRC(1)
    figX.plotPartialRCSupply(1)
    figX.plotPartialRCDemand(1)
    
    




