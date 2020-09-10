
typedef std::vector<TreeNode*> HaplotypeDescription;

class Tree
{
private:

	struct PQcompare
	{
		bool operator()(const std::pair<uint32_t, uint32_t>& a, const std::pair<uint32_t, uint32_t>& b)
		{
		    return a.first > b.first;
		}
	};

	// Arrays of nodes
	// lastGenerationHaplotypes, currentGenerationHaplotypes, ancestralGenerationHaplotypes and ancestral states all work as
	//     x[patch_index][ind_index * 2 + haplo_index]
	std::vector<std::vector<HaplotypeDescription>> lastGenerationHaplotypes;     
	std::vector<std::vector<HaplotypeDescription>> currentGenerationHaplotypes;

	// Security checks 
	uint32_t nbNodesInTree;
	bool currentStatesGotJustComputed;
	bool isCurrentlyBuildingHaplotype;

	// For optimization
	uint32_t TotalNbNodesCurrentGeneration;
	uint32_t TotalNbIndividualsCurrentGeneration;

	// While producing new haplotypes
	HaplotypeDescription newHaplotypeDescription;
	uint32_t parent_patch_index_newHaplotypeDescription;
	uint32_t parent_ind_index_newHaplotypeDescription;

	// outputing the tree
	OutputFile* file;
	
	


	void addNode(TreeNode* node);	
	std::queue<TreeNode*> initiateFIFO(std::vector<std::vector<HaplotypeDescription>>& oneGeneration, bool onlyNodesWithoutChildren);
	void pruneFromFIFO(std::queue<TreeNode*>& FIFO);
	std::pair<int,int> findNodeIDin(TreeNode*& node, std::vector<std::vector<HaplotypeDescription>>& listOfNodes);
	void destroyAllNodesFromGenerationAndParents(std::vector<std::vector<HaplotypeDescription>>& generationNodes);
	void swapLastAndcurrentGenerationHaplotypes();
	void destroyNode(TreeNode* node);
	void computeCurrentStates();
	void addMutationsToNode(TreeNode* node);
	void mergeCurrentHaplotypes();
	HaplotypeDescription mergeHaplotype(HaplotypeDescription& original);
	void destroyListOfNodes(std::vector<TreeNode*>& nodes);
	std::vector<TreeNode*> getRootsAndSetChildren();
	std::vector<TreeNode*> placeMutationsOnTree(std::queue<TreeNode*> FIFO, TreeNode* lastOfGenerationNode);
	std::string getNodeName(std::map<TreeNode*, std::string>& nodeNames, TreeNode* node, uint32_t& serialNumber);

	double computeMeanCoalescenceTime(std::vector<TreeNode*> offsprings);
	double computeTt();
	double computeTt_fast();
	double computeTs();
	double computeTs_fast();
	uint32_t computeCoalescentTimeBetweenTwoNodes(TreeNode* X, TreeNode* Y);
	double computeTb();
	double computeTb_fast();


public:

	void initialize();
	Tree();
	~Tree();
	
	void newGeneration();
	void addChildHaplotype_setParentInfo(uint32_t parent_patch, uint32_t parent_ind_index);
	void addChildHaplotype_addNode(uint32_t haploIndex, uint32_t from, uint32_t to);
	void addChildHaplotype_finished(uint32_t child_patch_index);


	void pruneDeadLineages();
	std::vector<std::vector<std::vector<std::vector<bool>>>> getCurrentStates();
	std::vector<std::vector<double>> getCurrentStates_frequencies();
	void indicateOutputFile(OutputFile* f);
	void printToFile(std::queue<TreeNode*> FIFO);

	std::vector<double> computeCoalescenceFstStatistics();
};

