/*

 Author: Remi Matthey-Doret

    MIT License

    Copyright (c) 2017 Remi Matthey-Doret

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.


Note of things to improve:
    Test how much slower it is to have N=1e5, L=10 vs N=10, L=1e5 to estimate the cost of having Individuals not contiguous in memory

    Improve performance of keywords (rep, seq, ...)

    Improve run function in Rwrapper for cases where the input is too big to be given through the command line (shell's limit)

    When using several environment, fitnessMap should take migration rate into account. This migration rate can vary through time and therefore fitnessMap should be redefined for faster simulations

 */

//#define CALLENTRANCEFUNCTIONS  

// Libraries
#include <iostream>
#include <list>
#include <stack>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <math.h>
#include <queue>
#include <ctime>
#include <set>
#include <array>
#include <chrono>
#include <thread>
//#include <boost/algorithm/string.hpp>
//#include <boost/tokenizer.hpp>
#include <cstdlib>
#include <stdlib.h>
//#include <boost/cstdint.hpp> // Fixed size integers for FDR
#include <signal.h>
#include <unordered_map>
#include "compressedString.cpp"
//#include <moreRNG/FastDiceRoller.cpp>
//#include <moreRNG/xorshiro128.cpp>

typedef float T3type;
typedef double fitnesstype;

// forward class declaration and global variables
#include "ForwardClassDeclarations.h"
#include "GlobalVariables.h"

// RNG
#include "RNG_wrapper.h"
#include "RNG_wrapper.cpp"


#include "T1_locusDescription.h"
#include "T1_locusDescription.cpp"

// variousSmallFunctions.cpp
#include "variousSmallFunctions.cpp"

// other .h
#include "BinaryFileToRead.h"
#include "CompressedSortedDeque.h" // also contains the .h for CompressedSortedDeque::CompressedSortedDequeBlock and CompressedSortedDeque::iterator 
#include "Genealogy.h"
#include "Option.h"
#include "OptionContainer.h"
#include "GeneralParameters.h"
#include "T56_memoryManager/T56_memoryManager.h"
#include "Walker.h"
#include "GeneticSampler.h"
//#include "SimulationTracker.h"
#include "DispersalData.h" 
#include "ResetGenetics.h"
#include "OutputFile.h"
#include "TreeNode.h"
#include "Tree.h"
#include "KillOnDemand.h"
#include "GeneticMap/GeneticMap.h"
#include "SampleSequenceDataContainer/SampleSequenceData.h"
#include "SampleSequenceDataContainer/SampleSequenceDataContainer.h"
#include "T4TreeRec/T4TreeRec.h"
#include "T8TreeRecording/MemoryPool.h"
#include "T8TreeRecording/T8Segment.h"
#include "T8TreeRecording/T8TreeRecording.h"

#include "ForcedMigration.h"
#include "T7stuff/T7Parameters.h"
#include "T7stuff/T7Gene.h"
#include "T7stuff/develop.h"
#include "SpeciesSpecificParameters.h"
#include "AllParameters.h"
#include "FromLocusToTXLocusElement.h"
#include "Pop.h"
#include "ZipIterator.h"
#include "Haplotype.h"
#include "Individual.h"
#include "InputReader.h"
#include "LifeCycle.h"
#include "LifeCycle_InnerDataClasses.h"
#include "OutputWriter.h"
#include "Patch.h"


// .cpp
#include "BinaryFileToRead.cpp"
#include "T7stuff/T7Gene.cpp"
#include "T7stuff/develop.cpp"
#include "CompressedSortedDeque.cpp"      // also contains the .cpp for CompressedSortedDeque::iterator and CompressedSortedDeque::CompressedSortedDequeBlock
#include "Genealogy.cpp"
#include "Option.cpp"
#include "OptionContainer.cpp"
#include "AllParameters.cpp"
#include "DispersalData.cpp"
#include "ResetGenetics.cpp"
#include "FromLocusToTXLocusElement.cpp"
#include "GeneralParameters.cpp"
#include "T56_memoryManager/T56_memoryManager.cpp"
#include "Walker.cpp"
#include "GeneticSampler.cpp"
#include "Pop.cpp"
#include "ZipIterator.cpp"
#include "Haplotype.cpp"
#include "Individual.cpp"
#include "InputReader.cpp"
#include "LifeCycle_InnerDataClasses.cpp"
#include "LifeCycle.cpp"
#include "OutputFile.cpp"
#include "OutputWriter.cpp"
#include "Patch.cpp"
//#include "SimulationTracker.cpp"
#include "KillOnDemand.cpp"
#include "GeneticMap/GeneticMap.cpp"
#include "SampleSequenceDataContainer/SampleSequenceData.cpp"
#include "SampleSequenceDataContainer/SampleSequenceDataContainer.cpp"
#include "ForcedMigration.cpp"
#include "SpeciesSpecificParameters.cpp"

#include "T8TreeRecording/MemoryPool.cpp"
#include "T8TreeRecording/T8Segment.cpp"
#include "T8TreeRecording/T8TreeRecording.cpp"

#include "T4TreeRec/Segment.cpp"
#include "T4TreeRec/ReversibleVector.cpp"
#include "T4TreeRec/Edge.cpp"
#include "T4TreeRec/EdgeTable.cpp"
#include "T4TreeRec/HaplotypeOfSegments.cpp"
#include "T4TreeRec/HaplotypesContainer.cpp"
#include "T4TreeRec/LastGenerationIDmap.cpp"
#include "T4TreeRec/NodeGenetics.cpp"
#include "T4TreeRec/NodeTable.cpp"
#include "T4TreeRec/T4TreeRec.cpp"


#include "codeToInsert.cpp"
    
//#include <boost/tokenizer.hpp> // was used to split strings but now the code only uses STL functions
//#include <omp.h> // If cannot find this file, then the Alejandro answer at 'http://stackoverflow.com/questions/35134681/installing-openmp-on-mac-os-x-10-11' can be helpful!


/*
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # ############################################################# #
 # ############################################################# #
 # ###################### Global variables ##################### #
 # ############################################################# #
 # ############################################################# #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 */

// SimBit Version
std::string getSimBitVersionLogo()
{
    std::string VERSION("version 4.16.25");
    std::string s;
    s.reserve(250);
    s += "\t  ____  _           ____  _ _   \n";
    s += "\t / ___|(_)_ __ ___ | __ )(_) |_ \n";
    s += "\t \\___ \\| | '_ ` _ \\|  _ \\| | __|\n";
    s += "\t  ___) | | | | | | | |_) | | |_ \n";
    s += "\t |____/|_|_| |_| |_|____/|_|\\__| ";
    s += VERSION + "\n\n"; //+ "\tWARNING: This is a version that will need further debugging of new features before being released\n\n" ;
    return s;
}

// allSpecies
std::vector<std::pair<Pop,Pop>> allSpecies;

// Global allParameters
AllParameters allParameters;

// global pointers to attributes of allParameters
SpeciesSpecificParameters* SSP  = nullptr;
GeneralParameters* GP           = nullptr;   // Also contains patchSizes for all species (duplicate info than from SSP)

//OutputWriter
OutputWriter outputWriter;

#ifdef DEBUG
int NBMIGRANTS = 0;
int NBT1MUTATIONS = 0;
#endif


/*
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # ############################################################# #
 # ############################################################# #
 # ##################### Initialize statics #################### #
 # ############################################################# #
 # ############################################################# #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 */

std::vector<T3type> Individual::T3_IndPhenotype;
std::vector<std::vector<double>> Individual::T7phenotypeOverTime;
std::vector<double> Individual::T7_IndPhenotype;
std::vector<fitnesstype> Individual::fitnessComponents = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // This is necessary for when there is no selection at all.
std::string OutputFile::GeneralPath;
std::vector<int> Pop::T2LociToCorrect;
std::string OutputFile::sequencingErrorStringToAddToFilnames = "";
std::vector<uint32_t> LifeCycle::recombination_breakpoints = {INT_MAX};
LifeCycle::ParentsData LifeCycle::PD;
bool KillOnDemand::justAboutToKill = false;
std::vector<int> T4TreeRec::generationsToKeepInTheTree;
std::vector<double> T7Gene::K_values_map = {5, 36.94528, 272.99075, 2017.14397, std::numeric_limits<double>::max()};
const long double Pop::limitMinimumSumOfFitnessAuthorized = 1e-300;//std::numeric_limits<long double>::min()*262144;
std::vector<uint32_t> T8TreeRecording::staticGeneticData = {};

//size_t T8Segment::static_mutate_position;
//size_t T8Segment::static_mutate_newsize;
//size_t T8Segment::static_mutate_mut;


/*
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # ############################################################# #
 # ############################################################# #
 # ######################## MAIN funtion ####################### #
 # ############################################################# #
 # ############################################################# #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 */

int main(int argc, char *argv[])
{

/*
###########################################
############## INITALIZATION ##############
###########################################
*/

    std::cout << getSimBitVersionLogo();

#ifdef DEBUG
    std::cout << "Run in DEBUG mode" << std::endl;
#endif

#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Run in CALLENTRANCEFUNCTIONS mode" << std::endl;
#endif  

    // measure time of process
    outputWriter.SetBeforeInitializationTime();

    // Keep track of the number of times a T2 block neede to be corrected
    unsigned long long int nbT2LociCorrections = 0;

    // Set GP*
    GP  = allParameters.getGPaddress();

    // Set Parameters
    allParameters.SetParameters(argc, argv);

    if (outputWriter.LogfileType != 0)
    {
        std::cout << "\t\tLogfile path: " << outputWriter.get_OutputFiles(Logfile)[0].getPath() << "\n\n";
    }
        

#ifdef DEBUG
    std::cout << "-------- Parameters read -------- " << std::endl;
#endif


    // start at the desired generation
    for (int GenerationChangeIndex = 0 ; GenerationChangeIndex < GP->__GenerationChange.size() ; GenerationChangeIndex++ )
    {
        int generation = GP->__GenerationChange[GenerationChangeIndex];
        if (generation > GP->startAtGeneration)
        {
            break;
        } else
        {
            allParameters.UpdateParameters();
        }
    }

    GP->CurrentGeneration = GP->startAtGeneration; 
    assert(GP->nbSpecies > 0);
    assert(SSP == nullptr);

    /*
    GP->CurrentGeneration = GP->startAtGeneration;
    allParameters.UpdatePopsAndParameters();
    */

    // Build population and do stuff
    {
        RNG_type rngForReset = GP->rngw.getRNG(); // reset after building populations. This is to avoid potential problem when loading seed but one of the species need random numbers to be drawn to be build.
        allSpecies.reserve(GP->nbSpecies);
        for (int speciesIndex = 0; speciesIndex < GP->nbSpecies ; speciesIndex++)
        {
            #ifdef DEBUG
            std::cout << "-------- species " << speciesIndex << " --------" << std::endl;
            #endif

            // set SSP*
            SSP = allParameters.getSSPaddress(speciesIndex);
            
            
            // Create populations
            Pop pop_odd(SSP->readPopFromBinary);


            SSP->resetGenetics.resetPopIfNeeded(pop_odd);
            if (SSP->Gmap.T56_nbLoci)
            {
                pop_odd.toggleT56MutationsIfNeeded();
            }

            #ifdef DEBUG
                std::cout << "-------- pop has been created --------" << std::endl;
            #endif

            //SSP->simTracker.Initialization(pop_odd); // The initialization will differ depending on whether the starting population explicity has not variance or whether it does.

            //#ifdef DEBUG
            //    std::cout << "-------- simTracker initialized --------" << std::endl;
            //#endif

            pop_odd.CalculateFitnesses();

            #ifdef DEBUG
                std::cout << "-------- Fitness of pop computed --------" << std::endl;
            #endif

            /*
            Pop pop_even;
            pop_even = pop_odd;
            if (GP->startAtGeneration + 1 <= GP->nbGenerations) // If something is simulated (a user might want nothing to be simualated just to extract a binary file)
            {
                assert(pop_odd.CumSumFits.size() == GP->PatchNumber);
                assert(pop_odd.CumSumFits[0].size() > 0);
                assert(pop_even.CumSumFits.size() == GP->PatchNumber);
                assert(pop_even.CumSumFits[0].size() > 0);
            }
                

            #ifdef DEBUG
                std::cout << "-------- pop has been copied --------" << std::endl;
            #endif


            //SSP->ClearT1_Initial_AlleleFreqs(); // We don't need the vector 'T1_Initial_AlleleFreqs' Anymore
            allSpecies.push_back(std::pair<Pop,Pop>(std::move(pop_odd), std::move(pop_even)));
            */

            allSpecies.push_back(std::pair<Pop,Pop>(
                pop_odd,
                std::move(pop_odd)
            ));

            #ifdef DEBUG
                std::cout << "-------- pop even and odd added to 'allSpecies' --------" << std::endl;
            #endif

            GP->CurrentGeneration = GP->startAtGeneration; // Not needed for the moment but it is just in case I modify the ::initialize method
            if (SSP->Gmap.T4_nbLoci) SSP->T4Tree.initialize(allSpecies[speciesIndex].second);
            if (SSP->Gmap.T8_nbLoci) SSP->T8Tree.initialize(allSpecies[speciesIndex].second);

            #ifdef DEBUG
                std::cout << "-------- T4Tree initialized --------" << std::endl;
            #endif

            allSpecies[speciesIndex].second.PrintBinaryFile();    // Save population in binary file if you asked for it. It will also save the seed to binary if SSP points to the parameters for the last speciesIndex 
            if (GP->burnInUntilT4Coal_check_every_N_generations == -1) outputWriter.WriteOutputs(allSpecies[speciesIndex].second); // WriteOutputs if you asked for it

            #ifdef DEBUG
            std::cout << "-------- Headers for output created --------" << std::endl;
            #endif
        }
        GP->rngw.setRNG(rngForReset);
    }
        

    // Security as we are outside a species specific context
    SSP          = nullptr;

    // Stop here if it is a dryRun
    if (GP->DryRun)
    {
        std::cout << "You asked for a DryRun. SimBit will therefore abort\n";
        abort();
    }

    // Time calculation
    outputWriter.SetBeforeSimulationTime();
    outputWriter.PrintInitializationTimeToLogFile();
    



    /*
    #############################################
    ############## Generation loop ##############
    #############################################
    */
        

    // Loop through each generation
    std::vector<bool> isFirstPopOffspring(GP->nbSpecies, false);

    if (GP->burnInUntilT4Coal_check_every_N_generations == -1)
    {
        GP->CurrentGeneration = GP->startAtGeneration + 1;
    } else
    {
        GP->CurrentGeneration = std::numeric_limits<int>::lowest() + 1;
    }

    for ( ; GP->CurrentGeneration <= GP->nbGenerations ;  GP->CurrentGeneration++ )
    {

        #ifdef DEBUG
        std::cout << "NBMIGRANTS = " << NBMIGRANTS << " ";
        std::cout << "NBT1MUTATIONS = " << NBT1MUTATIONS << "\n";
        NBMIGRANTS = 0;
        NBT1MUTATIONS = 0;
        #endif
        outputWriter.PrintGeneration();

        // Change pops and parameters if needed. Will loop through each species and will use the globals variables.
        if (GP->CurrentGeneration > GP->startAtGeneration + 1) allParameters.UpdatePopsAndParameters(); // The 'if' is here in case we want a neutral burnIn        
            

        // Set patch size after the end of the current generation. Is used for species interactions
        GP->setAllPatchSizePreviousGenerationIfNeeded();

        // Loop through each species to do the breeding
        for (int speciesIndex = 0; speciesIndex < GP->nbSpecies ; speciesIndex++)
        {
            // Get the big pointer for the species parameters
            SSP = allParameters.getSSPaddress(speciesIndex);

            if (SSP->whenDidExtinctionOccur == -1) // is not exctinct yet
            {
                assert(SSP->nbSubGenerationsPerGeneration > 0);
                for (int subGeneration = 0 ; subGeneration < SSP->nbSubGenerationsPerGeneration ; subGeneration++ )
                {
                    // get Pops
                    isFirstPopOffspring[speciesIndex]     = !isFirstPopOffspring[speciesIndex];
                    Pop& pop_Offspring  = isFirstPopOffspring[speciesIndex] ? allSpecies[speciesIndex].first : allSpecies[speciesIndex].second;
                    Pop& pop_Parent     = isFirstPopOffspring[speciesIndex] ? allSpecies[speciesIndex].second : allSpecies[speciesIndex].first;                    

                    // Early Forced migration
                    SSP->forcedMigration.forceMigrationIfNeeded(pop_Parent, true);

                    LifeCycle::BREEDING_SELECTION_DISPERSAL(pop_Offspring, pop_Parent);

                    if (GP->CurrentGeneration == GP->nbGenerations && subGeneration == SSP->nbSubGenerationsPerGeneration-1)
                    {
                        //std::cout << "About to free parent population memory\n";
                        pop_Parent.freeMemory();
                        if (SSP->Gmap.T56_nbLoci) pop_Offspring.shrink_to_fitT56();
                        //std::cout << "Parent population memory has been freed\n";
                    }

                    // In case, the T2 blocks capacity have been overpassed, correct them or abort.
                    nbT2LociCorrections += pop_Offspring.correctT2Loci();
                }
            } 
            // Check the T2 blocks approximation is not too awefully slow
            if (nbT2LociCorrections > (50 + GP->nbSpecies * GP->CurrentGeneration / 2))
            {
                std::cout << "There are a lot of T2 blocks correction for reaching the maximal number of mutations (255 mutations per block per individual). There were " << nbT2LociCorrections << " corrections in only " << GP->CurrentGeneration << " generations (for a simulation with " << GP->nbSpecies << " species)! SimBit will abort.\n"; 
                abort();
            }   
        }
        // Security as we are outside a species specific context
        SSP          = nullptr;

        // Loop through each species to update patchSize, to resetGenetics and to write outputs
        std::vector<bool> neutralBurnIn_hasSpeciesCoalesced(GP->nbSpecies, false);
        for (int speciesIndex = 0; speciesIndex < GP->nbSpecies ; speciesIndex++)
        {
            // Get the big pointers for the species
            SSP = allParameters.getSSPaddress(speciesIndex);

            // get Pop
            Pop& pop_Offspring  = isFirstPopOffspring[speciesIndex] ? allSpecies[speciesIndex].first : allSpecies[speciesIndex].second;

            if (GP->CurrentGeneration >= 0) // Don't get there if in burnIn
            {
                // Redefine individual types
                SSP->redefineIndividualTypesIfNeeded(pop_Offspring);

                // Kill individuals
                if (SSP->fecundityForFitnessOfOne != -1)
                {
                    SSP->killIndividualsIfAskedFor();
                }

                // reset genetics if needed (under user requirement only)
                SSP->resetGenetics.resetPopIfNeeded(pop_Offspring);

                // Code inserted
                codeToInsert(pop_Offspring);

                // Late Forced migration
                SSP->forcedMigration.forceMigrationIfNeeded(pop_Offspring, false);

                // WriteOutputs if you asked for it
                outputWriter.WriteOutputs(pop_Offspring); 

                // Save population in binary file if you asked for it
                pop_Offspring.PrintBinaryFile();


                // killOnDemand
                if (SSP->killOnDemand.mustKill(pop_Offspring))
                {
                    outputWriter.WriteOutputs(pop_Offspring);
                    std::string forSpeciesText;
                    if (GP->nbSpecies > 1)
                    {
                        forSpeciesText = " for species " + SSP->speciesName;
                    }
                    std::cout << "\n\n\t\tSimulation has been killed on demand"<< forSpeciesText << " at generation "<< GP->CurrentGeneration <<"!\n\n";
                    goto endOfSimulation;
                }
            }


            // Burn in -- that code is a bit long and ugly for being in main.cpp. Sorry
            else 
            {
                GP->testIfEndOfT4BurnIn(pop_Offspring, neutralBurnIn_hasSpeciesCoalesced);
            }
        }

        // Security as we are outside a species specific context
        SSP          = nullptr;

        // Check the simulation is not doing something crazy stupidly slow with blcok T2
        if (nbT2LociCorrections > 0)
        {
            std::cout << "\n\t\tThere was a need to correct T2 blocks for an excess of mutations on the block. Each correction is quite a slow process so you should attempt to keep the number of correction to the a low value (better just never reach 255 mutations per individual per block). You might want to adjust your parameters. There were " << nbT2LociCorrections << " corrections in total";
            abort();
        }
    }

    

    /*
    ##########################################
    ############## Finalization ##############
    ##########################################
    */
       
    endOfSimulation: 
    //std::cout << "\n\n" << std::flush;

    // Time calculation
    outputWriter.SetAfterSimulationTime();
    outputWriter.PrintSimulationTimeToLogFile();
    

    // Print a 'Good Job' message on the Logfile.
    outputWriter.FinalizeLogfile();

    // Print a 'Good Job' message on the SO.
    std::cout << "\t\t   --> Simulation is over. Good Job! <--" << "\n" << std::endl;


    return 0;
}
