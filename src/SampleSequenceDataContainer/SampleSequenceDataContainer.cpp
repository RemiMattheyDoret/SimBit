void SampleSequenceDataContainer::start()
{
	SSD_index = 0;
}

SampleSequenceData& SampleSequenceDataContainer::getNextSequence()
{
	assert(SSD_index < SSDs.size());
	return SSDs[SSD_index++];
}

bool SampleSequenceDataContainer::areThereMoreSequences() const
{
	return SSD_index < SSDs.size();
}

void SampleSequenceDataContainer::readInput(InputReader& input)
{
	while (input.IsThereMoreToRead() && input.PeakNextElementString() != "time")
	{
		SSDs.push_back(input);
	}
}