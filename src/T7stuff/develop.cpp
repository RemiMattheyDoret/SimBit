


void Individual::prepareDevelop(std::vector<std::vector<OneProtEffect>>& protEffects, std::vector<OneProtEffect>& basicSignalEffects)
{
	// Linearize the genes in an array (instead of two as there is currently one per haplotype)
    /*
	std::vector<T7Gene&> genes;
	for (unsigned haploIndex = 0 ; haploIndex < 2 ; ++haploIndex )
	{
		Haplotype& haplo = getHaplo(haploIndex);
		for (unsigned geneInHaploIndex = 0 ; geneInHaploIndex < haplo.nbT7Genes() ; ++geneInHaploIndex)
		{
            auto& gene = haplo.getT7_Allele(geneInHaploIndex);
			genes.push_back(gene);
		}
	}
    */

    size_t nbGenes = getHaplo(0).nbT7Genes() + getHaplo(1).nbT7Genes();
    
    assert(protEffects.size() == nbGenes);
    assert(basicSignalEffects.size() == nbGenes);
	

	// Cis and trans effects
	for (unsigned bigRecipientGeneIndex = 0 ; bigRecipientGeneIndex < nbGenes ; ++bigRecipientGeneIndex)
	{
        // get recipient gene
        bool recipientHaploIndex = bigRecipientGeneIndex >= getHaplo(0).nbT7Genes();
        size_t recipientGeneIndex;
        if (recipientHaploIndex) recipientGeneIndex = bigRecipientGeneIndex; else recipientGeneIndex = bigRecipientGeneIndex - getHaplo(0).nbT7Genes();
        auto& recipientGene = getHaplo(recipientHaploIndex).getT7_Allele(recipientGeneIndex);
		for (unsigned bigCausalGeneIndex = 0 ; bigCausalGeneIndex < nbGenes ; ++bigCausalGeneIndex)
		{
            // get causal gene
            bool causalHaploIndex = bigCausalGeneIndex >= getHaplo(0).nbT7Genes();
            size_t causalGeneIndex;
            if (causalHaploIndex) causalGeneIndex = bigCausalGeneIndex; else causalGeneIndex = bigCausalGeneIndex - getHaplo(0).nbT7Genes();
            auto& causalGene    = getHaplo(causalHaploIndex).getT7_Allele(causalGeneIndex);


            double cisEffect = recipientGene.cis[causalGene.targetCisSite].cisEffect;
            auto k = recipientGene.cis[causalGene.targetCisSite].k;
            assert(k < T7Gene::K_values_map.size());
            double K_val = T7Gene::K_values_map[k];
            double transEffect = causalGene.trans;

			bool isEnhancer = cisEffect > 0.0;
			double effectMagnitude = 1 - exp(fabs(cisEffect * transEffect));
			assert(effectMagnitude > 0.0);

            if (K_val != std::numeric_limits<double>::max() && effectMagnitude != 0.0)
            {
                protEffects[bigRecipientGeneIndex].push_back({bigCausalGeneIndex, isEnhancer, effectMagnitude, K_val, causalGene.phenotypicEffects});
            }
		}

        {
            double cisEffect = recipientGene.cis[0].cisEffect;
            auto k = recipientGene.cis[0].k;
            assert(k < T7Gene::K_values_map.size());
            double K_val = T7Gene::K_values_map[k];
            double transEffect = SSP->T7devpars.basic_signal_trans_effect;

            bool isEnhancer = cisEffect > 0.0;
            double effectMagnitude = 1 - exp(fabs(cisEffect * transEffect));
            assert(effectMagnitude > 0.0);
            effectMagnitude *= SSP->T7devpars.basic_signal_conc; // Direclty multiply byu the concentration as it is a constant
            
            std::vector<PhenotypicEffect> pE;
            basicSignalEffects[bigRecipientGeneIndex] = {0, isEnhancer, effectMagnitude, K_val, pE};
        }
            
	}
}




void Individual::develop(const int& Habitat)
{
    (void) Habitat;
	////////////////////
	/// Preparations ///
	////////////////////


	// Get total number of genes
	unsigned nbgenesInGenome = haplo0.nbT7Genes() + haplo1.nbT7Genes();
	assert(nbgenesInGenome);

	// Prepare develop
	std::vector<size_t> RNAconc(nbgenesInGenome, 0);   // RNA concs for each gene in the genome. They represent number of elemntss and not really a concentration
	std::vector<size_t> protconc(nbgenesInGenome, 0);  // prot concs for each gene in the genome. They represent number of elemntss and not really a concentration
	std::vector<std::vector<OneProtEffect>> protEffects(nbgenesInGenome); // who from (index in genes) and effect (cis and trans combined) // protEffects[recipient] -> returns a vector of affecting genes that are coded with three values (unsigned, the index in genes of the affecting gene; bool, wheter or not is enhancing; double, the magnitude of the effect)
    std::vector<OneProtEffect> basicSignalEffects(nbgenesInGenome);
	prepareDevelop(protEffects, basicSignalEffects); // will fill up 'RNAconc' and 'protconc' by reference



	//////////////////////////////////
	/// Big while loop preparation ///
	//////////////////////////////////

	unsigned age = 0;
	unsigned deltaT;
	std::vector<double> transcriptionRates(nbgenesInGenome);
	std::vector<double> oldTranscriptionRates(nbgenesInGenome);

    // If there is a need to record how the phenotype would evolve over time
    T7phenotypeOverTime.resize( SSP->T7phenpars.nbDimensions);
    T7_IndPhenotype.resize( SSP->T7phenpars.nbDimensions);
    for (size_t dim = 0 ; dim < SSP->T7phenpars.nbDimensions ; ++dim)
    {
        T7phenotypeOverTime[dim].resize(SSP->T7phenpars.agesAtwhichPhenotypeIsSampled.size());
        T7_IndPhenotype[dim] = 0.0;
    }
    

    ///////////////////////
    /// Big while loop  ///
    ///////////////////////
    size_t ageSampledIndex = 0;
	while (age < SSP->T7devpars.maxAge)
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

            // Basic signal effect. Note: magnitude of basic signal is already ultiplied by concentration
            if (basicSignalEffects[geneIndex].K_val != std::numeric_limits<double>::max())
            {
                totalProtein += SSP->T7devpars.basic_signal_conc;
                if (basicSignalEffects[geneIndex].isEnhancer)
                {
                    enhancing *= 1 - basicSignalEffects[geneIndex].magnitude / (totalProtein + basicSignalEffects[geneIndex].K_val);
                } else
                {
                    repressing *= 1 - basicSignalEffects[geneIndex].magnitude / (totalProtein + basicSignalEffects[geneIndex].K_val);
                }
            }

            // Protein effects
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

            transcriptionRates[geneIndex] = SSP->T7devpars.basal_transcription_rate * (1 - enhancing) * repressing;
        }


        /*
        /////////////////////////////////////////
		/// Change abundances of basic signal ///
		/////////////////////////////////////////

        // decay
        if (SSP->T7devpars.stochasticDevelopment)
        {
            std::poisson_distribution<double> d(protein_decayRate * deltaT * basic_signal_conc);
            basic_signal_conc -= d(GP->mt);
        } else
        {
            basic_signal_conc -= protein_decayRate * deltaT * basic_signal_conc;
        }
    	
        // new basical_signal
        basic_signal_conc += SSP->T7devpars.basic_signal_conc_rateInput * deltaT; // could add stoachsticity in that too!

        // Make sure it is positive
        if(basic_signal_conc < 0) basic_signal_conc = 0;
        */

        /////////////////////////////////////////////////////////////////////////////////////
		/// Change abundances of mRNAs and proteins and change change in phenotypic trait ///
		/////////////////////////////////////////////////////////////////////////////////////
        
        for(int geneIndex = 0; geneIndex < nbgenesInGenome; ++geneIndex)
        {
            // RNA concentration
            {
                auto rate = deltaT *
                        (
                            transcriptionRates[geneIndex]            // input
                            -
                            RNAconc[geneIndex] * SSP->T7devpars.mRNA_decayRate // decay
                        );
                if (SSP->T7devpars.stochasticDevelopment)
                {
                    RNAconc[geneIndex] -= GP->rngw.poisson(rate);
                } else
                {
                    RNAconc[geneIndex] -= rate;
                }
            }

            // protein concentration
            {
                auto rate = deltaT * 
                    (
                        RNAconc[geneIndex] * SSP->T7devpars.translationRate // input
                        -
                        protconc[geneIndex] * SSP->T7devpars.protein_decayRate // decay
                    );
                if (SSP->T7devpars.stochasticDevelopment)
                {
                    protconc[geneIndex] -= GP->rngw.poisson(rate);
                } else
                {
                    protconc[geneIndex] -= rate;
                }
            }

            // phenotypic traits
            {
                for (auto& geneEffect : protEffects[geneIndex])
                {
                    for (auto& phenoEffect : geneEffect.phenotypicEffects)
                    {
                        this->T7_IndPhenotype[phenoEffect.phenotypeDimension] += phenoEffect.magnitude * protconc[geneIndex];
                    }
                }
            }
        }

        /////////////////////////////////////
        /// reset old transcription rates ///
        /////////////////////////////////////

        transcriptionRates.swap(oldTranscriptionRates);


        ////////////////////////
        /// Record phenotype ///
        ////////////////////////

        assert(ageSampledIndex < SSP->T7phenpars.agesAtwhichPhenotypeIsSampled.size());
        auto nextAgeSampled = SSP->T7phenpars.agesAtwhichPhenotypeIsSampled[ageSampledIndex];
        if (age == nextAgeSampled)
        {
            for (size_t dim = 0 ; dim < SSP->T7phenpars.nbDimensions ; ++dim)
            {
                T7phenotypeOverTime[ageSampledIndex][dim] = this->T7_IndPhenotype[dim];
            }
            ++ageSampledIndex;
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
            auto oldDeltaT = deltaT;
            deltaT = SSP->T7devpars.maxDeltaT; // minimize deltaT. deltaT = max_time. Would probably be way too big!
            for(unsigned geneIndex = 0; geneIndex < nbgenesInGenome; geneIndex++)
            {
                double diff = fabs(oldTranscriptionRates[geneIndex] - transcriptionRates[geneIndex]);
                if (diff > 0)
                {
                    unsigned deltaT_tmp = SSP->T7devpars.EPSILON * oldDeltaT * oldTranscriptionRates[geneIndex] / diff;
                    if (deltaT_tmp < deltaT)
                        deltaT = deltaT_tmp;
                }
            }

    
            // Check for slowdown condition (protein levels far from equilibrium)
            for(unsigned geneIndex = 0; geneIndex < nbgenesInGenome; geneIndex++)
            {
                if( fabs(RNAconc[geneIndex] * SSP->T7devpars.translationRate - protconc[geneIndex] * SSP->T7devpars.protein_decayRate) / (1 + protconc[geneIndex]) > 0.25)
                {
                    deltaT /= 2;
                    break;
                }
            }


            // Make sure deltaT is plausible. Jeremy used SSP->T7devpars.maxDeltaT=5. I used SSP->T7devpars.maxDeltaT=6. I have no justification for either of these values
            assert(deltaT > 0);
            if (deltaT >= SSP->T7devpars.maxDeltaT)
            {
                std::cout << "In Individual::develop, estimated a deltaT of " << deltaT << ", which is greater than 6. Sounds a bit too big. Something might going wrong!\n";
                abort();
            }
        }

        
        if (age + deltaT > nextAgeSampled)
        {
            deltaT = nextAgeSampled - age;
            age = nextAgeSampled;
        } else
        {
            age += deltaT;
        }
	}
    assert(age == SSP->T7devpars.maxAge);
}
