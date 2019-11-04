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

// No need to include as everything is compiled together.
// Do not rename this file as it is directly included from main.cpp
// You have access to two global pointers, SSP and GP, the pointer to the Species Specific Parameters and to the General Parameters, respectively
// 'codeToInsert' is run for every species and at every (sub)generation.
// 'codeToInsert' is run after resetGenetics
// 'codeToInsert' can typically be used to kill the job depending on the state of the population (such as wheather a given locus is fixed for example)

void codeToInsert(Pop& pop)
{
	/////////////////
	//// DEFAULT ////
	/////////////////

	/* By default there is not need to write anything in here*/

	
	/////////////////////////////////////////////////////
	//// Figure out what species we are dealing with ////
	/////////////////////////////////////////////////////

	// Use 'SSP->speciesName' (and 'SSP->speciesIndex')

	///////////////////////////////////////
	//// Run code some generation only ////
	///////////////////////////////////////
	
	/*
	auto generation = GP->CurrentGeneration;
	int runCodeEveryNGeneration = 10; // runs code every 10th generation
	if (generation % runCodeEveryNGeneration == 0)
	{
		// Write code here
	}
	*/

	///////////////////////////////////
	//// Compute alleles frequency ////
	///////////////////////////////////
	
	/*
	// For T5 loci, you can direclty use the methods:
			computeT56ntrlFrequencies
			computeT56selFrequencies
			computeT56Frequencies
			computePatchSpecificT56ntrlFrequencies
			computePatchSpecificT56selFrequencies
			computePatchSpecificT56Frequencies
		For example do
			auto T56alleleFrequencies = pop.computeT56Frequencies();
	
	// For other loci, you can compute it yourself. For example

	// Sum up number of '1' (or true) alleles
	std::vector<double> T1allFreqs(GP->PatchNumber, 0.0); // GP->PatchNumber is the number of patches
	for (unsigned patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		Patch& patch = pop.getPatch(patch_index);
		for (unsigned ind_index = 0 ; ind_index < SSP->patchSizes[patch_index] ; ++ind_index)
		{
			Individual& individual patch.getInd(ind_index);

			for (unsigned locus = 0 ; locus < SSP->T1_nbBits ; ++locus)
			{
				T1allFreqs[locus] += individual.haplo0.getT1_Allele(locus);
				T1allFreqs[locus] += individual.haplo1.getT1_Allele(locus);
			}	
		}
	}

	// Divide by the total number of possible alleles
	for (unsigned locus = 0 ; locus < SSP->T1_nbBits ; ++locus)
	{
		T1allFreqs[locus] /= 2.0 * SSP->totalPatchSize;
	}

	// Of course you can also compute the frequency of only the locus of interest


	*/
}

