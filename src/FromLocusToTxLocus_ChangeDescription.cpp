class FromLocusToTxLocus_ChangeDescription
{
	unsigned char locusThatWillChange;
	/*
		locusThatWillChange values meaning
			0 -> no meaning
			1 -> T1
			2 -> T2
			3 -> T3
			4 -> T4
			5 -> T56ntrl
			6 -> T56sel
	*/
	size_t T1;
	size_t T2;
	size_t T3;
	size_t T4;
	size_t T56ntrl;
	size_t T56sel;
};

class FromLocuToTxLocus
{
	std::vector<FromLocusToTxLocus_ChangeDescription> locusAtChanges;
	std::vector<size_t> T1;
	std::vector<size_t> T2;
	std::vector<size_t> T3;
	std::vector<size_t> T4;
	std::vector<size_t> T56ntrl;
	std::vector<size_t> T56sel;
	std::vector<size_t> fromLocusToLastChangeIndex;

	template<typename INT>
	getT1(INT locus)
	{
		auto CI = fromLocusToLastChangeIndex[locus];
		return locusAtChanges[CI].locusThatWillChange == 1 ? T1[CI] + locusAtChanges[CI].T1 : T1[CI];
	}

	template<typename INT>
	getT2(INT locus)
	{
		auto CI = fromLocusToLastChangeIndex[locus];
		return locusAtChanges[CI].locusThatWillChange == 2 ? T2[CI] + locusAtChanges[CI].T2 : T2[CI];
	}

	template<typename INT>
	getT3(INT locus)
	{
		auto CI = fromLocusToLastChangeIndex[locus];
		return locusAtChanges[CI].locusThatWillChange == 3 ? T3[CI] + locusAtChanges[CI].T3 : T3[CI];
	}

	template<typename INT>
	getT4(INT locus)
	{
		auto CI = fromLocusToLastChangeIndex[locus];
		return locusAtChanges[CI].locusThatWillChange == 4 ? T4[CI] + locusAtChanges[CI].T4 : T4[CI];
	}

	template<typename INT>
	getT56ntrl(INT locus)
	{
		auto CI = fromLocusToLastChangeIndex[locus];
		return locusAtChanges[CI].locusThatWillChange == 5 ? T56ntrl[CI] + locusAtChanges[CI].T56ntrl : T56ntrl[CI];
	}

	template<typename INT>
	getT56sel(INT locus)
	{
		auto CI = fromLocusToLastChangeIndex[locus];
		return locusAtChanges[CI].locusThatWillChange == 6 ? T56sel[CI] + locusAtChanges[CI].T56sel : T56sel[CI];
	}
};