#include <string>
#include <vector>

template<typename mut_type, position_type>
class popGenData
{	
public:
	std::vector<std::vector<std::vector<mut_type>>> oA_alleles;  // data[patch_index][haplo_index][mutation_index]
	std::vector<std::vector<std::vector<mut_type>>> oB_alleles;  // data[patch_index][mutation_index][haplo_index]
	std::vector<position_type> SNPpositions;
	std::vector<size_t> patchIndices;	
	size_t maxPatchIndex;


private:
	size_t readingVCF_current_haplo_index = 0;
	
	initializationWhenReadingVCF(std::vector<size_t>& pI);

	void readingVCF_addAllele(mut_type allele);
	void readingVCF_addNewSNP(position_type pos);
	void skipElementsInISS(size_t n, std::istringstream& iss);

public:
	
	static void popGenData<mut_type, position_type> readVCF(std::string& path, std::vector<size_t>& patchIndices);
};
