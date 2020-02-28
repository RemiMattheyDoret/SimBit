
###################
### Information ###
###################

# This Rwrapper is very simple and is only meant to ease the management and creation of simulations. Please have a look at the manual section "R wrapper" for more information

#####################
### Sourcing code ###
#####################

# Make sure to set the path to the WRrapper directory in the SimBit directory
source("Rwrapper/SimulationClass.R")
source("Rwrapper/ParameterGrid.R")


###########################
### Constant parameters ###
###########################

pathForCommands = "" # Just write an appropriate path for where to write the commands

patchNumber = "1"
patchCapacity = "unif 100"
nbGenerations = "1e4"
mutationRate = "unif 1e-6"



#####################################
### Varying parameters parameters ###
#####################################

# Example: Full factorial design varying the number of loci (nbLoci) and the strength of selection (selectionStrength)

varying_nbLoci = c(8, 16, 32, 64, 128)
varying_selectionStrength = c(0, 0.001, 0.01)
varying_dominanceCoefficient = c(1, 0.5, 0.2)

#################################
### Create grid of parameters ###
#################################

# Example: Full factorial design varying the number of loci (nbLoci) and the strength of selection (selectionStrength)

# When using ParameterGrid() variables in different lists will expand in full factorial design. Variables found in the same list do not expand

# In the below example, the selection strength and the dominance of coefficient are not full factorial as they are in the same list. Instead, for s=0 is always associated to h=1, s=0.001 is associated to h=0.5 and s=0.01 is associated to h=0.2

PGrid = ParameterGrid(
	list(nbLoci = varying_nbLoci),
	list(
		selectionStrength = varying_selectionStrength,
		dominanceCoefficient = varying_dominanceCoefficient
	),
	list(
		patchNumber = patchNumber,
		patchCapacity = patchCapacity,
		nbGenerations = nbGenerations,
		mutationRate = mutationRate
	),
	bigID = "MyProject" # This is optional 
)

# Setting bigID is option. It allows to set an identifier for each row. Otherwise you could do
# PGrid$BigID = "MyProject"
# PGrid$smallID = 1:nrow(PGrid)
# PGrid$ID = paste0(PGrid$BigID, ".", PGrid$smallID)


###########################################
### Loop through the grid of parameters ###
###########################################

for (row in 1:nrow(PGrid))
{
	# # # # # # # # # #
	### Gather data ###
	# # # # # # # # # #
	# The function 'GetParameterGridData' sets for each column of the parameter grid, a variable with teh same name than this column and taking the value of the row 'row' of this column
	GetParameterGridData(PGrid, row)



	# # # # # # # # # # # # #
	### Create simulation ###
	# # # # # # # # # # # # #

	### Initialization of the command
	simulation = Simulation$new(pathForCommands, bigID, smallID)


	### Set values
	simulation$set("PN", patchNumber)
	simulation$set("N", patchCapacity)
	simulation$set("T1_mu", mutationRate)
	simulation$set("nbGens", nbGenerations)
	simulation$set("L", "T1", nbLoci)
	simulation$set("T1_fit", "unif 1", 1-dominanceCoefficient*selectionStrength, 1-selectionStrength) # Do not forget that using the multfit assumption (which we do not here) would make your simulations faster!
	

	##### If you want to check the input you do
	simulation$catCommand()

	#### If you want to check where this input is written do
	simulation$getCommandPath()

	#### If you want to directly run the simulation from R, do 
	simulation$run("SimBit") # Here write the path to the executable (or just SimBit if you placed it in a location pointed by the PATH variable in bash)
}

