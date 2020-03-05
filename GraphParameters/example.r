# This example shows you how to use the R6 class ParameterFile. That being said, it is quite simple to just copy-paste the parameter you want from the logfile type 2 rather than using my ParameterFile interface if you prefer.
# Note that the graph will likely not be good looking. You will need to play around with the ggnet function to get exactly the graph you want. The purpose of this script is just to show you a small example of what you could do to graph parameters that are being printed on the logfile type 2


###########################
##### Instal packages #####
###########################

needed_packages <- c("network", "ggplot2", "GGally")
need_to_install_packages = setdiff(needed_packages, rownames(installed.packages()))
if (length(need_to_install_packages) > 0) {
	cat(paste("I will try to install the following package(s):", paste(need_to_install_packages, collapse=", ")))
  install.packages(need_to_install_packages)
}
for (needed_package in needed_packages)
{
	require(needed_package, character.only=TRUE)
}


###################################
##### Load log file of type 2 #####
###################################

paramFile = ParameterFile$new("/Users/remi/Documents/Biologie/programming/PopGenSimulator/SimBit/log.log")

##########################
##### Get parameters #####
##########################

generations = paramFile$readParameter("__GenerationChange")  # Contains a vector of the generations at which there were temporal changes (with the @G marker)
dispMatrices = paramFile$readDispersalMatrixOverTime()       # Contains a list of dispersal matrix. Each element of the list is for a different @G marker
patchCapacities = paramFile$readParameter("__patchCapacity") # Contains a list of vector of patch capacity. Each element of the list is for a different @G marker

######################################################################################
##### Loop over the generations at which there were changes for those parameters #####
######################################################################################

pdf("DemographyGraphs.pdf")
for (generation_index in 1:length(generations))
{
	dispMatrix    = dispMatrices[[generation_index]]
	patchCapacity = patchCapacities[[generation_index]]
	generation    = generations[generation_index]

	#################
	##### Graph #####
	#################

	net = network(dispMatrix, directed = TRUE)
	#net %v% "patchID" = paste("patch ", 1:nrow(dispMatrix))
	print(ggnet(net, node.alpha=0.5, label = 1:nrow(dispMatrix), weight = patchCapacity, arrow.size = 8, arrow.gap = 0) + ggtitle(paste("From generation", generation)))
}
dev.off()


