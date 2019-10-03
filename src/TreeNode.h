
////////////////
/// TreeNode ///
////////////////

class TreeNode
{
public:
	// Will be used at all time
	TreeNode* parent;
	size_t nbChildren;
	size_t left;   // included
	size_t right;  // excluded

	// will remain empty until we computeStates
	std::vector<TreeNode*> children;
	std::vector<bool> genotype;  

	TreeNode(TreeNode* p, size_t l, size_t r);
	TreeNode(size_t l, size_t r);
	TreeNode(TreeNode* n);

	void assignGenotypeOfFalses();
	void announceDeathToChildren();
};


