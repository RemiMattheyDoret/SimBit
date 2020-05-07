
////////////////
/// TreeNode ///
////////////////

TreeNode::TreeNode(TreeNode* p, size_t l, size_t r) // normal constructor
:parent(p), nbChildren(0), left(l), right(r)
{
	parent->nbChildren++;
	// parent->children.push_back(this); No need to create that for the moment
#ifdef DEBUG
	assert(r <= SSP->T4_nbLoci);
#endif
}

TreeNode::TreeNode(size_t l, size_t r) // no parent constructor
:parent(nullptr), nbChildren(0), left(l), right(r)
{
	// parent->children.push_back(this); No need to create that for the moment
#ifdef DEBUG
	assert(r <= SSP->T4_nbLoci);
#endif
}


TreeNode::TreeNode(TreeNode* n) // copy constructor
:parent(n->parent), nbChildren(n->nbChildren), left(n->left), right(n->right)
{
	std::cout << "In TreeNode::TreeNode: I should not be using that copy constructor!\n";
	abort();
}

void TreeNode::assignGenotypeOfFalses()
{
	genotype.resize(SSP->T4_nbLoci, false);
}

void TreeNode::announceDeathToChildren()
{
	for (auto& child : children)
		child->parent = nullptr;
}

