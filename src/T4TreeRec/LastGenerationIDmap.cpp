
template<typename INT1, typename INT2>
T4TreeRec::LastGenerationIDmap::LastGenerationIDmap(INT1 nbNodes, INT2 totalPatchSize)
{
	assert(nbNodes >= totalPatchSize*2);
	min = nbNodes - totalPatchSize*2;
	newIDsForOldIDs.resize(totalPatchSize*2);
}

void T4TreeRec::LastGenerationIDmap::add(const uint32_t oldID, const uint32_t newID)
{
	auto i = oldID - min;
	assert(i < newIDsForOldIDs.size());
	newIDsForOldIDs[i] = newID;
}

uint32_t T4TreeRec::LastGenerationIDmap::getNewID(const uint32_t oldID) const
{
	auto i = oldID - min;
	assert(i >= 0 && i < newIDsForOldIDs.size());
	return newIDsForOldIDs[i];
}


