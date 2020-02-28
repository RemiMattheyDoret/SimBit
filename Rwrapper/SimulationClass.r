
###########################
##### Instal packages #####
###########################

needed_packages <- c("R6", "stringdist")
need_to_install_packages = setdiff(needed_packages, rownames(installed.packages()))
if (length(need_to_install_packages) > 0) {
	cat(paste("I will try to install the following package(s):", paste(need_to_install_packages, collapse=", ")))
  install.packages(need_to_install_packages)  
}
for (needed_package in needed_packages)
{
	require(needed_package, character.only=TRUE)
}
	


###################
##### R6Class #####
###################



Simulation = R6Class(
    ### Class name
    "Simulation",

    ### Private members
    private = list(
        commandPath = "",
        ID = "",
        shared = {
            env <- new.env()
            env$optionNames <- c("seed", "random_seed", "printProgress", "nbGens", "nbGenerations", "startAtGeneration", "nbThreads", "S", "species", "nbSubGens", "nbSubGenerations", "PN", "PatchNumber", "N", "patchCapacity", "H", "Habitats", "T5_approximationForNtrl", "T56_approximationForNtrl", "T5_fit", "T5_FitnessEffects", "T5_compressData", "T5_compress", "L", "Loci", "ploidy", "fec", "fecundityForFitnessOfOne", "m", "DispMat", "DispWeightByFitness", "gameteDispersal", "InitialpatchSize", "cloningRate", "selfingRate", "matingSystem", "LogfileType", "sequencingErrorRate", "additiveEffectAmongLoci", "selectionOn", "T1_mu", "T1_MutationRate", "T1_epistasis", "T1_EpistaticFitnessEffects", "T2_mu", "T2_MutationRate", "T4_mu", "T4_MutationRate", "T4_maxAverageNbNodesPerHaplotype", "T5_mu", "T5_MutationRate", "T5_toggleMutsEveryNGeneration", "T5_freqThreshold", "T5_frequencyThresholdForFlippingMeaning", "T3_mu", "T3_MutationRate", "T3_pheno", "T3_PhenotypicEffects", "T3_fit", "T3_FitnessLandscape", "T3_DN", "T3_DevelopmentalNoise", "T1_fit", "T1_FitnessEffects", "T2_fit", "T2_FitnessEffects", "r", "RecombinationRate", "recRateOnMismatch", "FitnessMapInfo", "indTypes", "individualTypes", "resetGenetics", "indIni", "individualInitialization", "T1_ini", "T1_Initial_AlleleFreqs", "T5_ini", "T5_Initial_AlleleFreqs", "GP", "GeneralPath", "T1_vcf_file", "T1_VCF_file", "T1_LargeOutput_file", "T1_AlleleFreq_file", "Log", "Logfile", "Logfile_file", "T1_MeanLD_file", "T1_LongestRun_file", "T1_HybridIndex_file", "T1_ExpectiMinRec_file", "T2_LargeOutput_file", "SaveBinary_file", "T3_LargeOutput_file", "T3_MeanVar_file", "fitness_file", "fitnessSubsetLoci_file", "fitnessStats_file", "T1_FST_file", "T1_FST_info", "extraGeneticInfo_file", "patchSize_file", "extinction_file", "genealogy_file", "coalesce", "shouldGenealogyBeCoalesced", "T4_LargeOutput_file", "T4_vcf_file", "T4_VCF_file", "T4_SFS_file", "T1_SFS_file", "T4_printTree", "T4_coalescenceFst_file", "T5_vcf_file", "T5_VCF_file", "T5_SFS_file", "T5_AlleleFreq_file", "T5_LargeOutput_file", "outputSFSbinSize", "eco", "speciesEcologicalRelationships", "popGrowthModel", "stochasticGrowth", "swapInLifeCycle", "Overwrite", "readPopFromBinary", "DryRun", "centralT1LocusForExtraGeneticInfo", "killOnDemand", "geneticSampling_withWalker", "individualSampling_withWalker")
            env
        }
    ),
    active = list(
        optionNames = function(value)
        {
            if (missing(value))
            {
                private$shared$optionNames
            } else
            {
                private$shared$optionNames <- value
            }
        }
    ),



    ### Public members
    public = list(
        initialize = function(pathToWriteCommandsTo, bigID, smallID)
        {
            private$ID = paste0(bigID, ".", smallID)
            path = pathToWriteCommandsTo

            if (substring(path, nchar(path), nchar(path)) == "/")
            {
                path = substring(path, 0, nchar(path)-1)
            }

            if (path == "")
            {
            	path = bigID
            } else
            {
            	path = paste0(path, "/", bigID)
            }

            dir.create(path, showWarnings = FALSE)
            private$commandPath <- paste0(path, "/command_", private$ID, ".txt")

            if (file.exists(private$commandPath))
            {
                file.remove(private$commandPath)
            }
            file.create(private$commandPath)
            
        },

        ### Set option
        set = function(optionName, ...)
        {
            if (!any(optionName == private$shared$optionNames))
            {
                distances = stringdist(optionName, private$shared$optionNames)
                alternatives = optionNames[which(order(distances) < 3)]

                stop(
                    paste0(
                    "\n\tReceived option name ", optionName, " to method 'set' for simulation ID ", private$ID, ".\n\tDid you maybe mean ", alternatives[1], ", ", alternatives[2], ", or ", alternatives[3], "?\n"
                    )
                )
            }

            entries = list(...)
            for (entry_index in 1:length(entries))
            {
                if (any(class(entries[[entry_index]]) == "data.frame"))
                {
                    entries[[entry_index]] = as.matrix(entries[[entry_index]])   
                }
                entries[[entry_index]] = paste(entries[[entry_index]], collapse=" ")
            }
            entry = paste(entries, collapse=" ")

            write(paste0("--", optionName, " ", entry), private$commandPath, append=TRUE)
        },

        ### print options
        printAllOptions = function()
        {
            print(private$shared$optionNames)
        },

        ### Print current command
        catCommand = function()
        {
            cat( readChar(private$commandPath, file.info(private$commandPath)$size) )
        },

        ### runSimulation
        run = function(executable)
        {
            system(paste(executable, "f", private$commandPath, "a"))
        },

        ### Get command path
        getCommandPath = function()
        {
        	return ( private$commandPath )
        }
    )
)

#### Example
# simulation = Simulation$new("Test","1","/Users/remi/test")
# simulation$set("PN", 1)
# simulation$set("N", "unif 100")
# simulation$set("L", "T1 8")
# simulation$set("T1_mu", "unif 1e-7")
# simulation$set("nbGens", "100")
# simulation$catCommand()
# simulation$getCommandPath()
# simulation$run("SimBit")

