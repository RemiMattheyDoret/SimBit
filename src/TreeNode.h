
////////////////
/// TreeNode ///
////////////////

class TreeNode
{
public:
	// Will be used at all time
	TreeNode* parent;
	uint32_t nbChildren;
	uint32_t left;   // included
	uint32_t right;  // excluded

	// will remain empty until we computeStates
	std::vector<TreeNode*> children;
	std::vector<bool> genotype;  

	TreeNode(TreeNode* p, uint32_t l, uint32_t r);
	TreeNode(uint32_t l, uint32_t r);
	TreeNode(TreeNode* n);

	void assignGenotypeOfFalses();
	void announceDeathToChildren();
};


