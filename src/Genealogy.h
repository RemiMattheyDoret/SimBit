class Genealogy // Contained in outputWrite
{
	// Format OffID_motherID_fatherID
	// example: P2I212_P2I12_P3I403
private:
	std::string buffer;
	std::string P = "P";
	std::string I = "I";
	std::string underscore = "_";
	std::string whitespace = " ";
	int generationFrom;
	int generationTo;
	bool isItTime=false;
	bool willItEverBeTime=false;
	int coalesceGenealogyFrequency;

	void setupBothParentTables(std::vector<std::vector<bool>>& parentsForNext, std::vector<std::vector<bool>>& parents);
	void resetParentsToFalse(std::vector<std::vector<bool>>& Parents);
	void getFields(std::vector<std::string>& fields, std::string& token, std::string filePath);
	void setParents(std::vector<std::string>& fields, std::vector<std::vector<bool>>& parents);
	void getParentsIDsFromFile(OutputFile& file, std::vector<std::vector<bool>>& parents);
	bool isIndividualParent(std::string& description, std::vector<std::vector<bool>>& parents);
	void rewriteFile(OutputFile& file, int generation, std::vector<std::vector<bool>>& Parents, std::vector<std::vector<bool>>& parentsForNext);
	OutputFile& modifyFilePath(OutputFile& file, int generation);
	
public:
	void printBufferIfIsTime();
	void startNewGeneration();
	void addOffspringIfIsTime(int offPatch, int offIndex, int motherPatch, int motherIndex, int fatherPatch, int fatherIndex);
	void coalesce(OutputFile& file);
	void mergeFiles(OutputFile& file);
	bool isTime();
	bool isTimeToCoalesce();
	void setGenealogyToNothing();
	void setGenealogyTimes(std::vector<int> times);
	void setcoalesceGenealogyFrequency(int x);
	bool isCoalesce();
	bool isTimeToMerge();
	void removeFilesAtStart(OutputFile& file);

	void writeOutWhoIsParent(std::vector<std::vector<bool>>& parents); // for debug purposes
};
