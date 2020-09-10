
template<typename mut_type, position_type>
void popGenData<mut_type, position_type>::initializationWhenReadingVCF(std::vector<size_t>& pI)
{
	patchIndices = pI;

	maxPatchIndex = 0;
	for (auto& patch_index : patchIndices)
		if (maxPatchIndex < patch_index) maxPatchIndex = patch_index;

	oA_alleles.resize(maxPatchIndex);
	oB_alleles.resize(maxPatchIndex);

	for (size_t patch_index = 0 ; patch_index < maxPatchIndex ; ++patch_index)
	{
		oA_alleles[patch_index].resize(patchIndices.size());
	}
}



	
template<typename mut_type, position_type>
void popGenData<mut_type, position_type>::readingVCF_addAllele(mut_type allele)
{
	assert(patchIndices.size() > readingVCF_current_haplo_index);
	auto patch_index = patchIndices[readingVCF_current_haplo_index];
	assert(data.oA_alleles[patch_index].size() > readingVCF_current_haplo_index);

	data.oA_alleles[patch_index][readingVCF_current_haplo_index].push_back(allele);
	data.oB_alleles[patch_index][SNPpositions.size()-1].push_back(allele);
	++readingVCF_current_haplo_index;
}

template<typename mut_type, position_type>
void popGenData<mut_type, position_type>::readingVCF_addNewSNP(position_type pos)
{
	SNPpositions.push_back(pos);
	data.oB_alleles.push_back({});
	readingVCF_current_haplo_index = 0;
}

