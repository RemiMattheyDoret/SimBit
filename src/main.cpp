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


Note for Remi of things to do:
    Test how much slower it is to have N=1e5, L=10 vs N=10, L=1e5 to estimate the cost of having Individuals not contiguous in memory

    When using several environment, fitnessMap should take migration rate into account. This migration rate can vary through time and therefore fitnessMap should be redefined for faster simulations

    As a general rule. Remi, please always remember that the interlocus n est l'interlocus qui vient juste apres le locus n. The last recombination position is therefore totalNbLoci - 2 as the last locus index is totalNbLoci - 1

 */

//#define CALLENTRANCEFUNCTIONS   

// Libraries
#include <iostream>
#include <list>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <sstream>
#include <iomanip>
#include <limits.h>
#include <algorithm>
#include <map>
#include <queue>
#include <ctime>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

/*#include <iostream>
#include <list>
#include <fstream>
#include <getopt.h>
#include <random>
#include <vector>
#include <string>
#include <cstring>
#include <limits>
#include <assert.h>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <math.h>
#include <limits.h>
#include <algorithm>
#include <map>
#include <unistd.h>
#include <ctime>
*/

// forward class declaration and global variables
#include "ForwardClassDeclarations.h"
#include "GlobalVariables.h"

// SmallOrderingSortingFunctions.cpp
#include "SmallOrderingSortingFunctions.cpp"

// other .h
#include "Genealogy.h"
#include "Option.h"
#include "OptionContainer.h"
#include "GeneralParameters.h"
#include "SimulationTracker.h"
#include "DispersalData.h" 
#include "ResetGenetics.h"
#include "TreeNode.h"
#include "Tree.h"
#include "SpeciesSpecificParameters.h"
#include "AllParameters.h"
#include "FromLocusToTXLocusElement.h"
#include "Pop.h"    
#include "Haplotype.h"
#include "Individual.h"
#include "InputReader.h"
#include "LifeCycle.h"
#include "OutputFile.h"
#include "OutputWriter.h"
#include "Patch.h"
#include "T1_locusDescription.h"

// .cpp    
#include "Genealogy.cpp"
#include "Option.cpp"
#include "OptionContainer.cpp"
#include "AllParameters.cpp"
#include "DispersalData.cpp"
#include "ResetGenetics.cpp"
#include "FromLocusToTXLocusElement.cpp"
#include "GeneralParameters.cpp"
#include "Pop.cpp"
#include "Haplotype.cpp"
#include "Individual.cpp"
#include "InputReader.cpp"
#include "LifeCycle.cpp"
#include "OutputFile.cpp"
#include "OutputWriter.cpp"
#include "Patch.cpp"
#include "SimulationTracker.cpp"
#include "SpeciesSpecificParameters.cpp"
#include "T1_locusDescription.cpp"
#include "TreeNode.cpp"
#include "Tree.cpp"
    
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
    std::string VERSION("version 4.0.5");
    std::string s;
    s.reserve(250);
    s += "\t  ____  _           ____  _ _   \n";
    s += "\t / ___|(_)_ __ ___ | __ )(_) |_ \n";
    s += "\t \\___ \\| | '_ ` _ \\|  _ \\| | __|\n";
    s += "\t  ___) | | | | | | | |_) | | |_ \n";
    s += "\t |____/|_|_| |_| |_|____/|_|\\__| ";
    s += VERSION + "\n\n";
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

std::vector<double> Individual::T3_IndPhenotype;
std::string OutputFile::GeneralPath;
std::vector<int> Pop::T2LociToCorrect;
std::string OutputFile::sequencingErrorStringToAddToFilnames = "";

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
        std::cout << "\t\tLogfile path:\n\t\t\t" << outputWriter.get_OutputFiles(Logfile)[0].getPath() << "\n\n";
    }
        

#ifdef DEBUG
    std::cout << "-------- Parameters read -------- " << std::endl;
#endif


    // start at the desired generation
    GP->UpdateParametersPatchNumber(0);
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


    // Build population and do stuff
    {
        std::mt19937 mtForReset = GP->mt; // reset after building populations. This is to avoid potential problem when loading seed but one of the species need random numbers to be drawn to be build.
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

            #ifdef DEBUG
                std::cout << "-------- pop has been created --------" << std::endl;
            #endif

            SSP->simTracker.Initialization(pop_odd); // The initialization will differ depending on whether the starting population explicity has not variance or whether it does.

            #ifdef DEBUG
                std::cout << "-------- simTracker initialized --------" << std::endl;
            #endif

            pop_odd.CalculateFitnesses();

            #ifdef DEBUG
                std::cout << "-------- Fitness of pop computed --------" << std::endl;
            #endif

            Pop pop_even = pop_odd;

            #ifdef DEBUG
                std::cout << "-------- pop has been copied --------" << std::endl;
            #endif

            #ifdef DEBUG
                std::cout << "-------- Initialization of simTracker terminated --------" << std::endl;
            #endif

            pop_odd.PrintBinaryFile();    // Save population in binary file if you asked for it. It will also save the seed to binary if SSP points to the parameters for the last speciesIndex 
            outputWriter.WriteOutputs(pop_odd); // WriteOutputs if you asked for it

            #ifdef DEBUG
            std::cout << "-------- Headers for output created --------" << std::endl;
            #endif


            //SSP->ClearT1_Initial_AlleleFreqs(); // We don't need the vector 'T1_Initial_AlleleFreqs' Anymore

            std::pair<Pop,Pop> popPair(std::move(pop_odd), std::move(pop_even));
            allSpecies.push_back(std::move(popPair));

            #ifdef DEBUG
                std::cout << "-------- pop even and odd added to 'allSpecies' --------" << std::endl;
            #endif
        }
        GP->mt = mtForReset;
    }
        

    // Security as we are outside a species specific context
    SSP          = nullptr;

    // Stop here if it is a dryRun
    if (GP->DryRun)
    {
        std::cout << "You asked for a DryRun.\n";
        abort();
    }

    // Time calculation
    outputWriter.SetBeforeSimulationTime();
    outputWriter.PrintInitializationTimeToLogFile();


    allParameters.UpdatePopsAndParameters();



    /*
    #############################################
    ############## Generation loop ##############
    #############################################
    */

    // Loop through each generation
    bool isGenerationOdd = false;
    for ( GP->CurrentGeneration = GP->startAtGeneration + 1 ; GP->CurrentGeneration <= GP->nbGenerations ;  GP->CurrentGeneration++ )
    {
        #ifdef DEBUG
        std::cout << "NBMIGRANTS = " << NBMIGRANTS << " ";
        std::cout << "NBT1MUTATIONS = " << NBT1MUTATIONS << "\n";
        NBMIGRANTS = 0;
        NBT1MUTATIONS = 0;
        #endif
        outputWriter.PrintGeneration();

        // Change pops and parameters if needed. Will loop through each species and will use the globals variables.
        allParameters.UpdatePopsAndParameters();

        // Loop through each species to do the breeding
        for (int speciesIndex = 0; speciesIndex < GP->nbSpecies ; speciesIndex++)
        {
            // Get the big pointer for the species parameters
            SSP = allParameters.getSSPaddress(speciesIndex);

            if (!SSP->isSpeciesExtinct)
            {
                assert(SSP->nbSubGenerationsPerGeneration > 0);
                for (int subGeneration = 0 ; subGeneration < SSP->nbSubGenerationsPerGeneration ; subGeneration++ )
                {
                    // get Pops
                    isGenerationOdd     = !isGenerationOdd;
                    Pop& pop_Offspring  = isGenerationOdd ? allSpecies[speciesIndex].first : allSpecies[speciesIndex].second;
                    Pop& pop_Parent     = isGenerationOdd ? allSpecies[speciesIndex].second : allSpecies[speciesIndex].first;

                    // Do selection dispersal, selection and breeding
                    LifeCycle::BREEDING_SELECTION_DISPERSAL(pop_Offspring,pop_Parent);

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
        for (int speciesIndex = 0; speciesIndex < GP->nbSpecies ; speciesIndex++)
        {
            // Get the big pointers for the species
            SSP = allParameters.getSSPaddress(speciesIndex);

            // get Pop
            Pop& pop_Offspring  = isGenerationOdd ? allSpecies[speciesIndex].first : allSpecies[speciesIndex].second;

            // reset genetics if needed (under user requirement only)
            SSP->resetGenetics.resetPopIfNeeded(pop_Offspring);

            // Save patchSize to GP
            GP->saveSSPPatchSize_toGP();

            // WriteOutputs if you asked for it
            outputWriter.WriteOutputs(pop_Offspring); // includes genealogy updates for coalescence

            // Save population in binary file if you asked for it
            pop_Offspring.PrintBinaryFile();

        }

        // Security as we are outside a species specific context
        SSP          = nullptr;

        // Check the simulation is not doing something crazy stupidly slow with blcok T2
        if (nbT2LociCorrections > 0)
        {
            std::cout << "\nThere was a need to correct T2 blocks for an excess of mutations on the block. Each correction is quite a slow process so you should attempt to keep the number of correction to the a low value (better just never reach 255 mutations per individual per block). You might want to adjust your parameters. There were " << nbT2LociCorrections << " corrections in total";
            abort();
        }
    }




    /*
    ##########################################
    ############## Finalization ##############
    ##########################################
    */
        
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
