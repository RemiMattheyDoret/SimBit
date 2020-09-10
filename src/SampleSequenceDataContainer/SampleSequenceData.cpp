void SampleSequenceData::start()
{
	if (lociList.size())
	{
		nextL = 0;
	} else
	{
		assert(to >= from);
		nextL = from;
	}
}

size_t SampleSequenceData::nbLoci() const
{
	if (lociList.size())
	{
		return lociList.size();
	} else
	{
		assert(to >= from);
		return to - from;
	}
}

bool SampleSequenceData::isMoreLoci()
{
	if (lociList.size())
	{
		return nextL < lociList.size();
	} else
	{
		return nextL < to;
	}
}

uint32_t SampleSequenceData::nextLocus()
{
	if (lociList.size())
	{
		return lociList[nextL++];
	} else
	{
		return nextL++;
	}
}


SampleSequenceData::SampleSequenceData(InputReader& input)
{
#ifdef DEBUG
	std::cout << "enters SampleSequenceData::SampleSequenceData(InputReader& input) with input = " << input.print() << "\n";
#endif


	// assert it starts with 'sequence'
	{
		std::string s = input.PeakNextElementString();
		if (s != "sequence")
		{
			std::cout << "In option --sampleSeq_file, expected the key 'sequence' but got '" << s << "' instead\n";
			abort();
		}
		input.skipElement();
	}
		

	patch = input.GetNextElementInt();
	if (patch >= GP->maxEverPatchNumber)
	{
		std::cout << "In option --sampleSeq_file, received a patch number (" << patch << ") that is greater than the highest number of patches ever in the simulation (" << GP->maxEverPatchNumber << ").\n";
		abort();
	}

	ind = input.GetNextElementInt();
	if (ind >= SSP->maxEverpatchCapacity[patch])
	{
		std::cout << "In option --sampleSeq_file, received a ind number (" << ind << ") that is greater than the highest capacity of this patch (patch "<<patch<<") ever in the simulation (" << SSP->maxEverpatchCapacity[patch] << ").\n";
		abort();
	}	

	haplo = input.GetNextElementInt();
	if (haplo != 0 && haplo != 1)
	{
		std::cout << "In option --sampleSeq_file, received a haplo number (" << haplo << ") that is not 0 or 1\n";
		abort();
	}	

	auto typeOfLociIndication = input.GetNextElementString();
	if (typeOfLociIndication == "fromto")
	{
		from = input.GetNextElementInt();
		to = input.GetNextElementInt();
		if (from > to || to >= SSP->Gmap.TotalNbLoci)
		{
			std::cout << "In option --sampleSeq_file, received the list of loci with the keyword 'fromto'. 'from' received is " << from << " and 'to' is " << to << ". Something is wrong. Either because one of these numbers is negative, because 'to' is actually not greater than 'from' or because one of these numbers is greater or equal to the total number of loci (which is equal to " << SSP->Gmap.TotalNbLoci << "\n";
			abort();
		}
	} else if (typeOfLociIndication == "list")
	{
		from = 0;
		to = 0;
		while (input.IsThereMoreToRead() && input.PeakNextElementString() != "sequence" && input.PeakNextElementString() != "time")
		{
			auto locus = input.GetNextElementInt();
			if (locus < 0 || locus > SSP->Gmap.TotalNbLoci)
			{
				std::cout << "In option --sampleSeq_file, received a locus (after keyword 'list') that is either lower than zero or greater to the total number of loci (TotalNbLoci = " << SSP->Gmap.TotalNbLoci << ").\n";
				abort();
			}
			lociList.push_back(locus);
		}
	} else
	{
		std::cout << "In option --sampleSeq_file, expected either keyword 'fromto' or 'list' but got "<<typeOfLociIndication<<" instead\n";
		abort();
	}
}
