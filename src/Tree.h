
typedef std::vector<TreeNode*> HaplotypeDescription;

class Tree
{
private:

	struct PQcompare
	{
		bool operator()(const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b)
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
	size_t nbNodesInTree;
	bool currentStatesGotJustComputed;
	bool isCurrentlyBuildingHaplotype;

	// For optimization
	size_t TotalNbNodesCurrentGeneration;
	size_t TotalNbIndividualsCurrentGeneration;

	// While producing new haplotypes
	HaplotypeDescription newHaplotypeDescription;
	size_t parent_patch_index_newHaplotypeDescription;
	size_t parent_ind_index_newHaplotypeDescription;
	
	


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

public:

	void initialize();
	Tree();
	~Tree();
	
	void newGeneration();
	void addChildHaplotype_setParentInfo(size_t parent_patch, size_t parent_ind_index);
	void addChildHaplotype_addNode(size_t haploIndex, size_t from, size_t to);
	void addChildHaplotype_finished(size_t child_patch_index);


	void pruneDeadLineages();
	std::vector<std::vector<std::vector<std::vector<bool>>>> getCurrentStates();
};

