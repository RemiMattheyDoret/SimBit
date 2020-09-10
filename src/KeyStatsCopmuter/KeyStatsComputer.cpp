void KeyStatsComputer::computeEHH_and_resetCloneList(PopGenData& data, std::deque<std::deque<size_t>>& cloneList, double& sumEHH, size_t& SNPi)
{
	// Reset all clone sets
	for (size_t cloneSet = 0 ; cloneSet < cloneList.size() ; )
	{
		// Split clone set
		std::deque<size_t> newCloneSet;
		for (size_t cloneIndex = 0 ; cloneIndex < cloneList[cloneSet].size() ; ++cloneIndex)
		{
			auto haplo = cloneList[cloneSet][cloneIndex];
			if (data.get(haplo, SNPi))
			{
				newCloneSet.push_back(haplo);
				cloneList[cloneSet].erase(cloneIndex);
			}
		}

		sumEHH += (double)(data.SNPposition(SNPi) - data.previousSNPposition(SNPi))/2 * pow((double) cloneList[cloneSet].size() / (double) data.nbHaplos(), 2) + pow((double) newCloneSet.size() / (double) data.nbHaplos(), 2);

		// Modify clone set
		if (cloneList[cloneSet].size() < 2 && newCloneSet.size() > 1)
		{
			// Switch clone set
			cloneList[cloneSet].swap(newCloneSet);
			++cloneSet
		} else
		{
			// Remove clone set if needed
			if (cloneList[cloneSet].size() < 2)
				cloneList.erase(cloneSet);
			else
				++cloneSet

			// Add clone set if needed
			if (newCloneSet.size() > 1)
				cloneList.push_back(newCloneSet);
		}
	}
}

std::pair<std::vector<double>, std::vector<double>> KeyStatsComputer::compute_iHS(PopGenData& data)
{
	std::vector<double> unstandard_iHS(data.nbSNPs());

	for (size_t SNPfocal = 0 ; SNPfocal < data.nbSNPs() ; ++SNPfocal)
	{
		double sumEHH_wt = 0.0;
		double sumEHH_mut = 0.0;

		std::deque<std::deque<size_t>> refCloneList_wt;
		std::deque<std::deque<size_t>> refCloneList_mut;
		for (size_t haplo = 0 ; haplo < data.nbHaplos() ; ++haplo)
		{
			if (data.get(haplo, SNPfocal))
				refCloneList_mut.push_back(haplo);
			else
				refCloneList_wt.push_back(haplo);
		}
		assert(refCloneList_wt.size());
		assert(refCloneList_mut.size());

		// left side
		if (SNPfocal > 1)
		{
			auto cloneLists_wt = refCloneList_wt;
			auto cloneLists_mut = refCloneList_mut;

			for (size_t SNPi = 0 ; SNPi < SNPfocal - 1 ; ++SNPi)
			{
				if (cloneLists_wt.size())computeEHH_and_resetCloneList(data, cloneLists_mut, sumEHH_mut, SNPi);
				if (cloneLists_mut.size())computeEHH_and_resetCloneList(data, cloneLists_wt, sumEHH_wt, SNPi);
				if (cloneLists_wt.size() == 0 && cloneLists_mut.size() == 0) break;
			}
		}

		// right side
		if (SNPfocal < data.nbHaplos()-1)
		{
			auto cloneLists_wt = refCloneList_wt;
			auto cloneLists_mut = refCloneList_mut;

			for (size_t SNPi = SNPfocal+1 ; SNPi < data.nbHaplos() ; ++SNPi)
			{
				if (cloneLists_wt.size())computeEHH_and_resetCloneList(data, cloneLists_mut, sumEHH_mut, SNPi);
				if (cloneLists_mut.size())computeEHH_and_resetCloneList(data, cloneLists_wt, sumEHH_wt, SNPi);
				if (cloneLists_wt.size() == 0 && cloneLists_mut.size() == 0) break;
			}
		}


		// Compute iHS
		unstandard_iHS[SNPfocal] = log(sumEHH_wt / sumEHH_mut);
	}

	// Standardize
	double mean_iHS = std::accumulate(unstandard_iHS.begin(), unstandard_iHS.end()) / unstandard_iHS.size();
	double var_iHS = 0.0;
	for (auto& uns_iHS : unstandard_iHS)
		var_iHS += pow(uns_iHS - mean_iHS, 2);
	var_iHS /= unstandard_iHS.size();

	std::vector<double> standard_iHS(data.nbSNPs());
	for (size_t SNPfocal = 0 ; SNPfocal < data.nbSNPs() ; ++SNPfocal)
		standard_iHS[SNPfocal] = (unstandard_iHS[SNPfocal] - mean_iHS) / var_iHS;

	return {unstandard_iHS, standard_iHS};
}


std::array<std::vector<double>, 3> KeyStatsComputer::compute_Hxy(PopGenData& data)
{
	for (size_t SNP = 0 ; SNP < data.nbSNPs() ; ++SNP)
	{
		double freq = 0.0;
		for (size_t haplo = 0 ; haplo < data.nbHaplos(); ++haplo)
		{
			freq += data.get(haplo, SNP);
		}
		freq /= data.nbHaplos();
		if (freq > 0.5)
		{
			H1 = freq * freq;
			H2 = (1-freq) * (1-freq);
		} else
		{
			H2 = freq * freq;
			H1 = (1-freq) * (1-freq);
		}
		H12 = ..

		Size of window??
		H1 = freq
	}
}


void KeyStatsComputer::computeAndWrite(PopGenData& data, std::string& outFilename)
{
	iHSData = compute_iHS(data);
	HxyData = compute_Hxy(data);
}



