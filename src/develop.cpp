To do:
	Needs to clean up the 'whoAffectsMe' of each gene every so often to remove gene IDs that dont exist anymore
	make separate functions the develop will call


class ContainerReturnedFromT7Develop
{
    std::vector<std::vector<double>> phenotypeOverTime; // a vector for each dimension
    std::vector<double> timesAtwhichPhenotypeIsSampled;
    ContainerReturnedFromT7Develop(std::vector<std::vector<double>> a, std::vector<double> b):phenotypeOverTime(a), timesAtwhichPhenotypeIsSampled(b) {}
};

class PhenotypicEffect
{
    unsigned phenotypeDimension;
    double magnitude;

    PhenotypicEffect(unsigned a, double b):phenotypeDimension(a), magnitude(b){}
};

class OneProtEffect
{
	unsigned indexCausal;
	bool isEnhancer;
	double magnitude;
	double K_val;
    std::vector<PhenotypicEffect>& phenotypicEffects;

	OneProtEffect(unsigned a, bool b, double c, char d, std::vector<PhenotypicEffect>& e):indexCausal(a), isEnhancer(b), magnitude(c), K_val(d), phenotypicEffects(e){}
};


void Individual::prepareDevelop(std::vector<double>& RNAconc, std::vector<double>& protconc, std::vector<std::vector<OneProtEffect>>& protEffects, const unsigned& nbgenesInGenome)
{
	// Linearize the genes in an array
	std::vector<T7Gene&> genes;
	for (unsigned haploIndex = 0 ; haploIndex < 2 ; ++haploIndex )
	{
		Haplotype& haplo = getHaplo(haploIndex);
		for (unsigned geneInHaploIndex = 0 ; geneInHaploIndex < haplo.T7Alleles.size() ; ++geneInHaploIndex)
		{
			genes.push_back(haplo.T7Alleles[geneInHaploIndex]);
		}
	}
	

	// Cis and trans effects
	for (unsigned recipientGeneIndex = 0 ; recipientGeneIndex < genes.size() ; ++recipientGeneIndex)
	{
		for (unsigned causalGeneIndex = 0 ; causalGeneIndex < genes.size() ; ++causalGeneIndex)
		{
			for (unsigned affectingGeneListIndex = 0 ; affectingGeneListIndex < genes[recipientGeneIndex].whoAffectsMe.size() ; ++affectingGeneListIndex)
			{
				if (genes[recipientGeneIndex].whoAffectsMe[affectingGeneListIndex] == genes[causalGeneIndex].ID)
				{
					// We found a gene that will affect recipient gene

					// Gather and compute effects
					double cisEffect = genes[recipientGeneIndex].howOthersAffectMe[affectingGeneListIndex].first;
					double K_val = valueForKmismatches[genes[recipientGeneIndex].howOthersAffectMe[affectingGeneListIndex].second];
					double transEffect = genes[causalGeneIndex].generalProtEffect;

					bool isEnhancer = cisEffect > 0.0;
					double effectMagnitude = 1 - exp(fabs(cisEffect * transEffect));
					assert(effectMagnitude > 0.0);

					// Make an entry in protEffects
					protEffects[recipientGeneIndex].push_back({causalGeneIndex, isEnhancer, effectMagnitude, K_val, genes[causalGeneIndex].phenotypicEffects})
				}
			}
		}
	}
	assert(genes.size() == protconc.size());
	assert(RNAconc.size() == protconc.size());
}



void Individual::develop(bool returnPhenotypeOverTime = false)
{
	////////////////////
	/// Preparations ///
	////////////////////

	// Hardcoded values
	const std::vector<double> valueForKmismatches = {5.00000, 36.94528, 272.99075, 2017.14397};


	// Get total number of genes
	unsigned nbgenesInGenome = haplo0.T7Alleles.size() + haplo1.T7Alleles.size();
	assert(nbgenesInGenome);

	// Prepare develop
	std::vector<double> RNAconc(nbgenesInGenome, 0.0);   // RNA concs for each gene in the genome
	std::vector<double> protconc(nbgenesInGenome, 0.0);  // prot concs for each gene in the genome
	std::vector<std::vector<OneProtEffect>> protEffects(nbgenesInGenome); // who from (index in genes) and effect (cis and trans combined) // protEffects[recipient] -> returns a vector of affecting genes that are coded with three values (unsigned, the index in genes of the affecting gene; bool, wheter or not is enhancing; double, the magnitude of the effect)
	prepareDevelop(RNAconc, protconc, protEffects, nbgenesInGenome); // will fill up 'RNAconc' and 'protconc' by reference



	//////////////////////////////////
	/// Big while loop preparation ///
	//////////////////////////////////

	unsigned age = 0;
	unsigned deltaT;
    unsigned SSP->maxDeltaT = 6;
	std::vector<double> transcriptionRates(nbgenesInGenome);
	std::vector<double> oldTranscriptionRates(nbgenesInGenome);

    // If there is a need to record how the phenotype would evolve over time
    std::vector<std::vector<double>> phenotypeOverTime; // a vector for each dimension
    std::vector<double> timesAtwhichPhenotypeIsSampled;
    if (SSP->T7_fitnessOverTime || returnPhenotypeOverTime)
    {
        phenotypeOverTime.reserve(SSP->T7_MaxAge / SSP->maxDeltaT * 1.2);
        timesAtwhichPhenotypeIsSampled.reserve(SSP->T7_MaxAge / SSP->maxDeltaT * 1.2);
    }

    ///////////////////////
    /// Big while loop  ///
    ///////////////////////

	while (age < SSP->T7_MaxAge)
	{
		//////////////////////////////////
		/// Compute transcriptionRates ///
		//////////////////////////////////

        for (unsigned geneIndex = 0 ; geneIndex < nbgenesInGenome ; ++geneIndex)
        {
        	// compute totalProtein
        	double totalProtein = 0.0;
        	for (auto& protEffect : protEffects[geneIndex])
        		totalProtein += protconc[protEffect.indexCausal];

        	// initalize enhancing and repressing
        	double enhancing = 1;
        	double repressing = 1;

        	for (auto& protEffect : protEffects[geneIndex])
        	{
        		if (protEffect.isEnhancer)
        		{
        			enhancing *= 1 - (protEffect.magnitude * protconc[protEffect.indexCausal]) / (totalProtein + protEffect.K_val);
        		} else
        		{
        			repressing *= 1 - (protEffect.magnitude * protconc[protEffect.indexCausal]) / (totalProtein + protEffect.K_val);
        		}
        	}
        	assert(repressing <= 1 && repressing >= 0);
            assert(enhancing <= 1 && enhancing >= 0);

            transcriptionRates[geneIndex] = SSP->basal_transcription_rate * (1 - enhancing) * repressing;
        }


        //////////////////
		/// Set detlaT ///
		//////////////////
		// detlaT represents the amoung by which 'age' will be advanced. The small deltaT is, the slower will the simulation be but the better the approximation.
        if (age == 0)
        {
        	deltaT = 1;
        } else
        {
        	// Find largest tolerable deltaT
        	deltaT = SSP->maxDeltaT; // minimize deltaT. deltaT = max_time. Would probably be way too big!
            for(unsigned geneIndex = 0; geneIndex < nbgenesInGenome; geneIndex++)
            {
            	double diff = fabs(oldTranscriptionRates[geneIndex] - transcriptionRates[geneIndex]);
            	unsigned deltaT_tmp;
                if (diff > 0)
                {
                	unsigned deltaT_tmp = SSP->T7_EPSILON * oldDeltaT * oldRates[geneIndex] / diff;
                	if (deltaT_tmp < deltaT)
                		deltaT = deltaT_tmp;
                }
            }

    
            // Check for slowdown condition (protein levels far from equilibrium)
            for(unsigned geneIndex = 0; geneIndex < nbgenesInGenome; geneIndex++)
            {
                if( fabs(RNAconc[geneIndex] * translation - protconc[geneIndex] * protein_decayRate) / (1 + protconc[geneIndex]) > 0.25)
                {
                    deltaT /= 2;
                    break;
                }
            }


            // Make sure deltaT is plausible. Jeremy used SSP->maxDeltaT=5. I used SSP->maxDeltaT=6. I have no justification for either of these values
            assert(deltaT > 0);
            if (deltaT >= SSP->maxDeltaT)
            {
            	std::cout << "In Individual::develop, estimated a deltaT of " << deltaT << ", which is greater than 6. Sounds a bit too big. Something might going wrong!\n";
            	abort();
            }
        }

        /////////////////////////////////////////
		/// Change abundances of basic signal ///
		/////////////////////////////////////////

        // decay
        if (SSP->T7_stochasticDevelopment)
        {
            std::poisson_distribution<double> d(protein_decayRate * deltaT * basic_signal_conc);
            basic_signal_conc -= d(GP->mt);
        } else
        {
            basic_signal_conc -= protein_decayRate * deltaT * basic_signal_conc;
        }
    	
        // new basical_signal
        basic_signal_conc += SSP->basic_signal_conc_rateInput * deltaT; // could add stoachsticity in that too!

        // Make sure it is positive
        if(basic_signal_conc < 0) basic_signal_conc = 0;


        /////////////////////////////////////////////////////////////////////////////////////
		/// Change abundances of mRNAs and proteins and change change in phenotypic trait ///
		/////////////////////////////////////////////////////////////////////////////////////
        
        for(int geneIndex = 0; geneIndex < nbgenesInGenome; ++geneIndex)
        {
            // RNA concentration
            {
                std::poisson_distribution<double> d(
                    deltaT *
                    (
                        transcriptionRates[geneIndex]            // input
                        -
                        RNAconc[geneIndex] * SSP->mRNA_decayRate // decay
                    )
                );
                RNAconc[geneIndex] -= d(GP->mt);
            }

            // protein concentration
            {
                std::poisson_distribution<double> d(
                    deltaT * 
                    (
                        RNAconc[geneIndex] * translationRate // input
                        -
                        protconc[geneIndex] * SSP->protein_decayRate // decay
                    )
                );
                protconc[geneIndex] -= d(GP->mt);
            }

            // phenotypic traits
            {
                for (auto& phenoEffect : protEffects[geneIndex].phenotypicEffects)
                {
                    this->T37_phenotypes[phenoEffect.phenotypeDimension] += phenoEffect.magnitude * protconc[geneIndex];
                }
            }
        }


        // If there is a need to record how the phenotype would evolve over time
        if (SSP->T7_fitnessOverTime || returnPhenotypeOverTime)
        {
            phenotypeOverTime.push_back(this->T37_phenotypes);
            timesAtwhichPhenotypeIsSampled.push_back(age);
        }
	}

    // add phenotype for last time point (whether or not other timepoints have been registered)
    assert(age == SSP->T7_MaxAge);
    phenotypeOverTime.push_back(this->T37_phenotypes);
    timesAtwhichPhenotypeIsSampled.push_back(age);
    assert(phenotypeOverTime.size() == timesAtwhichPhenotypeIsSampled.size());

    // return
    return {phenotypeOverTime, timesAtwhichPhenotypeIsSampled};

}
