GeneticSampler::GeneticSampler(){}



void GeneticSampler::set_nbRecombinations(const double& rate)
{
	recombinationtotalRate = rate;
	poisson_nbRecombinations =  std::poisson_distribution<size_t>(rate);
}

void GeneticSampler::set_T1_nbMuts(const double& rate)
{
	T1_mut_totalRate = rate;
	poisson_T1_nbMuts =  std::poisson_distribution<size_t>(rate);
}

void GeneticSampler::set_T2_nbMuts(const double& rate)
{
	T2_mut_totalRate = rate;
	poisson_T2_nbMuts =  std::poisson_distribution<size_t>(rate);
}

void GeneticSampler::set_T3_nbMuts(const double& rate)
{
	T3_mut_totalRate = rate;
	poisson_T3_nbMuts =  std::poisson_distribution<size_t>(rate);
}

void GeneticSampler::set_T56_nbMuts(const double& rate)
{
	T56_mut_totalRate = rate;
	poisson_T56_nbMuts =  std::poisson_distribution<size_t>(rate);
}


void GeneticSampler::set_recombinationPosition(const std::vector<double>& cumSumRates)
{
	assert(cumSumRates.size() == 1 || cumSumRates.size() == SSP->TotalNbLoci-1);
	if (SSP->TotalNbLoci > 1) 
	{
		/*
		std::cout << "recombinationtotalRate = " << recombinationtotalRate << "\n";
		std::cout << "cumSumRates.size() = " << cumSumRates.size() << "\n";
		std::cout << "cumSumRates.back() = " << cumSumRates.back() << "\n";
		*/
		isCstRate_recombinationPosition = false;
		if (SSP->TotalNbLoci <= 2)
		{
			// nothing to do
		} else if (cumSumRates.size() == 1) // meaning constant rate
		{
			isCstRate_recombinationPosition = true;
		} else 
		{
			assert(recombinationtotalRate == cumSumRates.back());
			if (SSP->geneticSampling_withWalker)
			{
				walker_recombinationPosition = Walker(cumSumRates);
			}
		}
		cumSumProbs_recombinationPosition = cumSumRates;
	}
}

void GeneticSampler::set_T1_mutationPosition(const std::vector<double>& cumSumRates)
{
	assert(cumSumRates.size() == 1 || cumSumRates.size() == SSP->T1_nbLoci);
	if (SSP->T1_nbLoci > 1)
	{
		isCstRate_T1_mutationPosition = false;
		if (SSP->T1_nbLoci == 1)
		{
			// nothing to do
		} else if (cumSumRates.size() == 1) // meaning constant rate
		{
			isCstRate_T1_mutationPosition = true;
		} else 
		{
			assert(T1_mut_totalRate == cumSumRates.back());
			if (SSP->geneticSampling_withWalker)
			{
				walker_T1_mutationPosition = Walker(cumSumRates);
			}
		}
		cumSumProbs_T1_mutationPosition = cumSumRates;
	}
}

void GeneticSampler::set_T2_mutationPosition(const std::vector<double>& cumSumRates)
{
	assert(cumSumRates.size() == 1 || cumSumRates.size() == SSP->T2_nbLoci);
	if (SSP->T2_nbLoci > 1)
	{
		isCstRate_T2_mutationPosition = false;
		if (SSP->T2_nbLoci == 1)
		{
			// nothing to do
		} else if (cumSumRates.size() == 1) // meaning constant rate
		{
			isCstRate_T2_mutationPosition = true;
		} else 
		{
			assert(T2_mut_totalRate == cumSumRates.back());
			if (SSP->geneticSampling_withWalker)
			{
				walker_T2_mutationPosition = Walker(cumSumRates);
			} 
		}
		cumSumProbs_T2_mutationPosition = cumSumRates;
	}
}

void GeneticSampler::set_T3_mutationPosition(const std::vector<double>& cumSumRates)
{
	assert(cumSumRates.size() == 1 || cumSumRates.size() == SSP->T3_nbLoci);
	if (SSP->T3_nbLoci > 1) 
	{
		isCstRate_T3_mutationPosition = false;
		if (SSP->T3_nbLoci == 1)
		{
			// nothing to do
		} else if (cumSumRates.size() == 1) // meaning constant rate
		{
			isCstRate_T3_mutationPosition = true;
		} else 
		{
			assert(T3_mut_totalRate == cumSumRates.back());
			if (SSP->geneticSampling_withWalker)
			{
				walker_T3_mutationPosition = Walker(cumSumRates);
			}
		}
		cumSumProbs_T3_mutationPosition = cumSumRates;
	}
}

void GeneticSampler::set_T56_mutationPosition(const std::vector<double>& cumSumRates)
{
	assert(cumSumRates.size() == 1 || cumSumRates.size() == SSP->T56_nbLoci);
	if (SSP->T56_nbLoci > 1)
	{
		isCstRate_T56_mutationPosition = false;
		if (SSP->T56_nbLoci == 1)
		{
			// nothing to do
		} else if (cumSumRates.size() == 1) // meaning constant rate
		{
			isCstRate_T56_mutationPosition = true;
		} else
		{
			assert(T56_mut_totalRate == cumSumRates.back());
			if (SSP->geneticSampling_withWalker)
			{
				walker_T56_mutationPosition = Walker(cumSumRates);
			}
		}
		cumSumProbs_T56_mutationPosition = cumSumRates;
	}
}





size_t GeneticSampler::get_nbRecombinations()
{
	if (recombinationtotalRate > 0.0)
	{
		return poisson_nbRecombinations(GP->rngw.getRNG());
	} else
	{
		return 0;
	}
		
}

size_t GeneticSampler::get_T1_nbMuts()
{
	if (T1_mut_totalRate > 0.0)
	{
		return poisson_T1_nbMuts(GP->rngw.getRNG());
	} else
	{
		return 0;
	}
		
}

size_t GeneticSampler::get_T2_nbMuts()
{
	if (T2_mut_totalRate > 0.0)
	{
		return poisson_T2_nbMuts(GP->rngw.getRNG());
	} else
	{
		return 0;
	}
		
}

size_t GeneticSampler::get_T3_nbMuts()
{
	if (T3_mut_totalRate > 0.0)
	{
		return poisson_T3_nbMuts(GP->rngw.getRNG());
	} else
	{
		return 0;
	}
		
}

size_t GeneticSampler::get_T56_nbMuts()
{
	if (T56_mut_totalRate > 0.0)
	{
		return poisson_T56_nbMuts(GP->rngw.getRNG());
	} else
	{
		return 0;
	}
		
}



size_t GeneticSampler::get_recombinationPosition()
{
	if (SSP->TotalNbLoci <= 2)
	{
		assert(SSP->TotalNbLoci==2);
		return 0;
	}
	if (isCstRate_recombinationPosition)
	{
		return GP->rngw.uniform_int_distribution(SSP->TotalNbLoci - 2);
	} else
	{
		if (SSP->geneticSampling_withWalker)
		{
			return walker_recombinationPosition(GP->rngw);
		} else
		{
			double rnd = GP->rngw.uniform_real_distribution(recombinationtotalRate);
			auto index = std::upper_bound(cumSumProbs_recombinationPosition.begin(), cumSumProbs_recombinationPosition.end(), rnd) - cumSumProbs_recombinationPosition.begin();
			assert(index >= 0 && index < SSP->TotalNbLoci - 1);
			return index;
		}
	}

}

size_t GeneticSampler::get_T1_mutationPosition()
{
	if (SSP->T1_nbLoci == 1)
	{
		return 0;
	}
	if (isCstRate_T1_mutationPosition)
	{
		return GP->rngw.uniform_int_distribution(SSP->T1_nbLoci - 1);
	} else
	{
		if (SSP->geneticSampling_withWalker)
		{
			return walker_T1_mutationPosition(GP->rngw);
		} else
		{
			double rnd = GP->rngw.uniform_real_distribution(T1_mut_totalRate);
			auto index = std::upper_bound(cumSumProbs_T1_mutationPosition.begin(), cumSumProbs_T1_mutationPosition.end(), rnd) - cumSumProbs_T1_mutationPosition.begin();
			assert(index >= 0 && index < SSP->T1_nbLoci);
			return index;
		}
	}

}

size_t GeneticSampler::get_T2_mutationPosition()
{
	if (SSP->T2_nbLoci == 1)
	{
		return 0;
	}
	if (isCstRate_T2_mutationPosition)
	{
		return GP->rngw.uniform_int_distribution(SSP->T2_nbLoci - 1);
	} else
	{
		if (SSP->geneticSampling_withWalker)
		{
			return walker_T2_mutationPosition(GP->rngw);
		} else
		{
			double rnd = GP->rngw.uniform_real_distribution(T2_mut_totalRate);
			auto index = std::upper_bound(cumSumProbs_T2_mutationPosition.begin(), cumSumProbs_T2_mutationPosition.end(), rnd) - cumSumProbs_T2_mutationPosition.begin();
			assert(index >= 0 && index < SSP->T2_nbLoci);
			return index;
		}
	}

}

size_t GeneticSampler::get_T3_mutationPosition()
{
	if (SSP->T3_nbLoci == 1)
	{
		return 0;
	}
	if (isCstRate_T3_mutationPosition)
	{
		return GP->rngw.uniform_int_distribution(SSP->T3_nbLoci - 1);
	} else
	{
		if (SSP->geneticSampling_withWalker)
		{
			return walker_T3_mutationPosition(GP->rngw);
		} else
		{
			double rnd = GP->rngw.uniform_real_distribution(T3_mut_totalRate);
			auto index = std::upper_bound(cumSumProbs_T3_mutationPosition.begin(), cumSumProbs_T3_mutationPosition.end(), rnd) - cumSumProbs_T3_mutationPosition.begin();
			assert(index >= 0 && index < SSP->T3_nbLoci);
			return index;
		}
	}

}

size_t GeneticSampler::get_T56_mutationPosition()
{
	if (SSP->T56_nbLoci == 1)
	{
		return 0;
	}
	if (isCstRate_T56_mutationPosition)
	{
		return GP->rngw.uniform_int_distribution(SSP->T56_nbLoci - 1);
	} else
	{
		if (SSP->geneticSampling_withWalker)
		{
			return walker_T56_mutationPosition(GP->rngw);
		} else
		{
			double rnd = GP->rngw.uniform_real_distribution(T56_mut_totalRate);
			auto index = std::upper_bound(cumSumProbs_T56_mutationPosition.begin(), cumSumProbs_T56_mutationPosition.end(), rnd) - cumSumProbs_T56_mutationPosition.begin();
			assert(index >= 0 && index < SSP->T56_nbLoci);
			return index;
		}
	}

}

