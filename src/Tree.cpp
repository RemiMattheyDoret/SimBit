
double Tree::computeMeanCoalescenceTime(std::vector<TreeNode*> offsprings)
{
	// Initialize variable to return
	double T = 0.0;

	// Security
	//uint32_t totalNbCoalescence = 0;
	double totalNbCoalescence = offsprings.size() - 1;

	// Initialize parents
	std::vector<TreeNode*> parents;
	parents.reserve(offsprings.size()); // Will never need to grow more

	// Loop through generations
	for (uint32_t generation = 1 ; true ; ++generation) // looping back in time
	{
		assert(generation <= GP->CurrentGeneration+1);

		// List parents
		for (uint32_t offspring_haplo_index = 0 ; offspring_haplo_index < offsprings.size() ; ++offspring_haplo_index)
		{
			parents.push_back(offsprings[offspring_haplo_index]->parent);
		}

		// Remove duplicate
		std::sort( parents.begin(), parents.end() );
		parents.erase( std::unique( parents.begin(), parents.end() ), parents.end() );


		// Add to T
		uint32_t nbCoalescence = offsprings.size() - parents.size();
		//totalNbCoalescence += nbCoalescence;
		/*
		std::cout << " ----------- \n";
		std::cout << "generation = " << generation << "\n";
		std::cout << "offsprings.size() = " << offsprings.size() << "\n";
		std::cout << "parents.size() = " << parents.size() << "\n";
		std::cout << "nbCoalescence = " << nbCoalescence << "\n";
		std::cout << "totalNbCoalescence = " << totalNbCoalescence << "\n";
		*/
		
		T += (double)(nbCoalescence * generation);

		// Condition for breaking and security
		assert(parents.size() <= offsprings.size());
		if (parents.size() == 1)
		{
			//std::cout << "Breaks after " << generation << " generations\n";
			break;
		} else
		{
			assert(parents.size() > 1);
		}

		// Clear offsprings and set parents as offsprings (parents become offsprings for next generation)
		offsprings.clear();
		parents.swap(offsprings);
	}

	/*
	std::cout << "parents.size() = " << parents.size() << "\n";
	std::cout << "offsprings.size() = " << offsprings.size() << "\n";
	std::cout << "totalNbCoalescence = " << totalNbCoalescence << "\n";
	std::cout << "expectedTotalNbCoalescence = " << expectedTotalNbCoalescence << "\n";
	*/
	//assert(totalNbCoalescence == expectedTotalNbCoalescence);

	return T / totalNbCoalescence;
}

double Tree::computeTt_fast()
{
	// Initialize offsprings
	std::vector<TreeNode*> offsprings;
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		uint32_t currentNbHaplos = currentGenerationHaplotypes[patch_index].size();
		offsprings.reserve(offsprings.size() + currentNbHaplos);

		for (uint32_t haplo_index = 0 ; haplo_index < currentNbHaplos ; ++haplo_index)
		{
			offsprings.push_back(currentGenerationHaplotypes[patch_index][haplo_index][0]);
		}
	}

	// computeMeanCoalescenceTime and cast it to double
	return computeMeanCoalescenceTime(offsprings);
}

double Tree::computeTs_fast()
{
	double Ts = 0.0;

	// Loop through patches
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		// Initialize offsprings
		std::vector<TreeNode*> offsprings;

		// Reserve
		uint32_t currentNbHaplos = currentGenerationHaplotypes[patch_index].size();
		offsprings.reserve(currentNbHaplos); // It will never need to grow

		// Gather offsprings
		for (uint32_t haplo_index = 0 ; haplo_index < currentNbHaplos ; ++haplo_index)
		{
			offsprings.push_back(currentGenerationHaplotypes[patch_index][haplo_index][0]);
		}


		// computeMeanCoalescenceTime
		Ts += computeMeanCoalescenceTime(offsprings);
	}

	// Cast to double, average, and return. This weight each patch equally. That might not be the best system!
	return Ts / (double) GP->PatchNumber;
}

double Tree::computeTs()
{
	double Ts = 0.0;
	for (uint32_t patch_i = 0 ; patch_i < GP->PatchNumber  ; ++patch_i)
	{		
		uint32_t Ts_patchi = 0;
		uint32_t nbComparisons = 0;
		for (uint32_t haplo_i = 0 ; haplo_i < currentGenerationHaplotypes[patch_i].size(); ++haplo_i)
		{
			for (uint32_t haplo_j = haplo_i ; haplo_j < currentGenerationHaplotypes[patch_i].size(); ++haplo_j)
			{
				++nbComparisons; // Even when 'patch_i == patch_j && haplo_i == haplo_j' it counts as a comparison for which coalescence time is zero
				Ts_patchi += computeCoalescentTimeBetweenTwoNodes
				(
					currentGenerationHaplotypes[patch_i][haplo_i][0],
					currentGenerationHaplotypes[patch_i][haplo_j][0]
				);
			}	
		}
		Ts += (double) Ts_patchi / (double) nbComparisons;
	}

	return Ts / (double) GP->PatchNumber;
}

uint32_t Tree::computeCoalescentTimeBetweenTwoNodes(TreeNode* X, TreeNode* Y)
{
	for (uint32_t generation = 0 ; true ; ++generation)
	{
		if (X == Y)
		{
			return generation;
		}

		X = X->parent;
		Y = Y->parent;
	}
}

double Tree::computeTb()
{
	uint32_t Tb = 0;
	uint32_t nbComparisons = 0;
	for (uint32_t patch_i = 0 ; patch_i < (GP->PatchNumber - 1)  ; ++patch_i)
	{
		for (uint32_t patch_j = patch_i + 1 ; patch_j < GP->PatchNumber ; ++patch_j)
		{
			for (uint32_t haplo_i = 0 ; haplo_i < currentGenerationHaplotypes[patch_i].size(); ++haplo_i)
			{
				for (uint32_t haplo_j = 0 ; haplo_j < currentGenerationHaplotypes[patch_j].size(); ++haplo_j)
				{
					++nbComparisons; // Even when 'patch_i == patch_j && haplo_i == haplo_j' it counts as a comparison for which coalescence time is zero
					Tb += computeCoalescentTimeBetweenTwoNodes
					(
						currentGenerationHaplotypes[patch_i][haplo_i][0],
						currentGenerationHaplotypes[patch_j][haplo_j][0]
					);
				}	
			}
		}
	}

	return (double) Tb / (double) nbComparisons;
}

double Tree::computeTb_fast()
{
	std::cout << "double Tree::computeTb_fast() is still undefined\n";
	abort();
}

double Tree::computeTt()
{
	uint32_t Tt = 0;
	uint32_t nbComparisons = 0;
	for (uint32_t patch_i = 0 ; patch_i < (GP->PatchNumber - 1)  ; ++patch_i)
	{
		for (uint32_t patch_j = patch_i ; patch_j < GP->PatchNumber ; ++patch_j)
		{
			for (uint32_t haplo_i = 0 ; haplo_i < currentGenerationHaplotypes[patch_i].size(); ++haplo_i)
			{
				for (uint32_t haplo_j = 0 ; haplo_j < currentGenerationHaplotypes[patch_j].size(); ++haplo_j)
				{
					++nbComparisons; // Even when 'patch_i == patch_j && haplo_i == haplo_j' it counts as a comparison for which coalescence time is zero
					Tt += computeCoalescentTimeBetweenTwoNodes
					(
						currentGenerationHaplotypes[patch_i][haplo_i][0],
						currentGenerationHaplotypes[patch_j][haplo_j][0]
					);
				}	
			}
		}
	}

	return (double) Tt / (double) nbComparisons;	
}


std::vector<double> Tree::computeCoalescenceFstStatistics()
{
	assert(SSP->TotalRecombinationRate == 0.0); // This assumption could be relaxed

	if (GP->PatchNumber > 1)
	{
		//pruneDeadLineages();

		double Ts = computeTs();
		std::cout << "Ts = " << Ts << std::endl;
		std::cout << "Ts_fast = " << computeTs_fast() << std::endl;
		double Tt = computeTt();
		std::cout << "Tt = " << Tt << std::endl;
		std::cout << "Tt_fast = " << computeTt_fast() << std::endl; 
		double Tb = computeTb();
		std::cout << "Tb = " << Tb << std::endl;
		std::cout << "Fst = " << (Tt - Ts) / Tt << std::endl;
		std::cout << "WCst = " << (Tb - Ts) / Tb << std::endl;
		//assert(Tt == Ttsecurity);

		return {Tt, Ts, Tb, (Tt - Ts) / Tt, (Tb - Ts) / Tb};
	} else
	{
		return {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
	}
}



std::string Tree::getNodeName(std::map<TreeNode*, std::string>& nodeNames, TreeNode* node, uint32_t& serialNumber)
{
	auto it = nodeNames.find(node);
	if (it == nodeNames.end())
	{
		// Not in map yet
		std::string r = std::to_string(serialNumber);
		nodeNames.insert(std::pair<TreeNode*,std::string>(node, r));
		serialNumber++;
		return r;
	} else
	{
		// Already in map
		return it->second;
	}
}


void Tree::printToFile(std::queue<TreeNode*> FIFO)
{
	std::map<TreeNode*, std::string> nodeNames;	
	uint32_t serialNumber = 0;

	std::string s;
	s += "Generation " + std::to_string(GP->CurrentGeneration) + "\n\n" + "Tree:\n";

	while (FIFO.size())
	{
		TreeNode* node = FIFO.front(); FIFO.pop();
		//assert(node->nbChilren == node->children.size());

		auto name = getNodeName(nodeNames, node, serialNumber);
		s += name + ":";
		for (TreeNode*& child : node->children)
		{
			s += " " + getNodeName(nodeNames, child, serialNumber);
		}
		s += "\n";
	}

	file->open();
	file->write(s);
	s.resize(0);

	s += "\ngenotypicRange:\n";
	for (auto it = nodeNames.begin() ; it != nodeNames.end() ; it++)
	{
		s += getNodeName(nodeNames, it->first, serialNumber) + ":" + std::to_string(it->first->left) + "-" + std::to_string(it->first->right) + "\n";
	}
	s += "\n";

	file->write(s);
	file->close();
}






void Tree::addNode(TreeNode* node)
{
	//std::cout << "\tadding new node from " << node->left << " to " << node->right << "\n";
	// The parent had its count of nbChildren increased directly by the constructor of TreeNode.
	newHaplotypeDescription.push_back(node);
	nbNodesInTree++;
}

void Tree::addChildHaplotype_finished(uint32_t child_patch_index)
{
	assert(currentGenerationHaplotypes.size() > child_patch_index);

#ifdef DEBUG
	{
		uint32_t expectedLeft = 0;

		for (uint32_t chunkIndex = 0 ; chunkIndex < newHaplotypeDescription.size() ; chunkIndex++)
		{
			std::cout << "Adding newHaplotypeDescription. Chunck " << chunkIndex << ": left = "<<newHaplotypeDescription[chunkIndex]->left<<", right = "<<newHaplotypeDescription[chunkIndex]->right<<"\n";
			assert(expectedLeft == newHaplotypeDescription[chunkIndex]->left);
			assert(expectedLeft <= newHaplotypeDescription[chunkIndex]->right);
			expectedLeft = newHaplotypeDescription[chunkIndex]->right;
		}
		std::cout << "expectedLeft = " << expectedLeft << "\n";
		assert(expectedLeft == SSP->Gmap.T4_nbLoci);
	}
#endif

	currentGenerationHaplotypes[child_patch_index].push_back(newHaplotypeDescription);
	TotalNbNodesCurrentGeneration += newHaplotypeDescription.size();
	TotalNbIndividualsCurrentGeneration++;
	newHaplotypeDescription.resize(0);

	assert(isCurrentlyBuildingHaplotype);
	isCurrentlyBuildingHaplotype               = false;
	//std::cout << "newHaplo added\n";
#ifdef DEBUG
	std::cout << "Exiting Tree::addChildHaplotype_finished\n";
#endif	
}

void Tree::pruneFromFIFO(std::queue<TreeNode*>& FIFO)
{
#ifdef DEBUG
	std::cout << "Enters pruneFromFIFO\n";
#endif
	while (FIFO.size() != 0)
	{
		TreeNode* node = FIFO.front(); 
		FIFO.pop();
		
		if (node->nbChildren != 0)
			std::cout << "node->nbChildren = " << node->nbChildren << "\n";
		assert(node->nbChildren == 0);

		TreeNode* parent  = node->parent;

		destroyNode(node); // destroy node before checking nbChildren of parent

		if (parent != nullptr) // if not an ancestral node
		{
			// parent->nbChildren--; This is done in destroyNode
		
			if (parent->nbChildren == 0)
			{
				FIFO.push(parent);
			}
		}
	}
}

void Tree::destroyNode(TreeNode* node)
{
	if (node->parent != nullptr)
	{
		assert(node->parent->nbChildren > 0);
		node->parent->nbChildren--;
	}
	
	nbNodesInTree--;
	delete node;
}

void Tree::initialize()
{
#ifdef DEBUG
  std::cout << "Enters Tree::initialize\n";
#endif		
	nbNodesInTree = 0;
	currentStatesGotJustComputed = true;
	isCurrentlyBuildingHaplotype = false;
	TotalNbIndividualsCurrentGeneration = 0;
	TotalNbNodesCurrentGeneration = 0;

	// Create first (parentless) nodes and put it in currentGenerationHaplotypes
	lastGenerationHaplotypes.resize(GP->PatchNumber);
	currentGenerationHaplotypes.resize(GP->PatchNumber);
	for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
	{
		for (int indHaplo_index = 0 ; indHaplo_index < (SSP->patchSize[patch_index] * 2) ; indHaplo_index++)
		{
			TreeNode* newNode = new TreeNode(0, SSP->Gmap.T4_nbLoci); // no parent constructor
			newNode->assignGenotypeOfFalses();
			// this->addChildHaplotype_setParentInfo(); no need to set parent info. Just set isCurrentlyBuildingHaplotype to true
			this->isCurrentlyBuildingHaplotype = true;
			this->addNode(newNode);
			this->addChildHaplotype_finished(patch_index);
			assert(!this->isCurrentlyBuildingHaplotype); // assert it is back to false
		}
	}


#ifdef DEBUG
  std::cout << "In Tree::initialize, "<<nbNodesInTree<<" nodes have been initializedTree()\n";
#endif	
}

void Tree::destroyListOfNodes(std::vector<TreeNode*>& nodes)
{
	for (auto& node : nodes)
		destroyNode(node);
}

void Tree::destroyAllNodesFromGenerationAndParents(std::vector<std::vector<HaplotypeDescription>>& generationNodes)
{
	auto FIFO = initiateFIFO(generationNodes, false); // false -> even nodes with children (nbChildren set to zero for those who have)
	pruneFromFIFO(FIFO);
}

Tree::Tree()
:nbNodesInTree(0), file(nullptr)
{}

void Tree::indicateOutputFile(OutputFile* f)
{
	file = f;
}

Tree::~Tree()
{
	if (file != nullptr)
	{
		if (!currentStatesGotJustComputed)
		{
			// initialize FIFO
			auto roots = getRootsAndSetChildren();
			std::queue<TreeNode*> FIFO;
			for (auto& root : roots)
			{	
				assert(root->nbChildren == root->children.size());
				FIFO.push(root);
			}

			printToFile(FIFO);
		}
	}

	if (nbNodesInTree > 0)
	{
		// destroy the whole tree // That's a bit of a funny way to destroy the tree. The current generation has no children and so everyone should be killed
		destroyAllNodesFromGenerationAndParents(currentGenerationHaplotypes);	
	}
	

	// check no nodes are left
	if (nbNodesInTree != 0)
	{
		std::cout << "Internal error in the destructor of Tree for T4. After attempt to remove all the nodes, it seems that there are " << nbNodesInTree << " nodes left, while it was supposed to have 0 nodes left!\n";
		abort();
	}

	// ugly pointer usage here
	delete file;
}

std::queue<TreeNode*> Tree::initiateFIFO(std::vector<std::vector<HaplotypeDescription>>& oneGeneration, bool onlyNodesWithoutChildren)
{
	std::queue<TreeNode*> FIFO;
	//std::cout << "oneGeneration.size() = " << oneGeneration.size() << "\n";
	assert(oneGeneration.size() == GP->PatchNumber);
	for (int patch_index = 0 ; patch_index < GP->PatchNumber; patch_index++)
	{
		for (auto& haploDescription : oneGeneration[patch_index])
		{
			for (auto& node : haploDescription)
			{
				if (onlyNodesWithoutChildren)
				{
					if (node->nbChildren == 0)
					{
						FIFO.push(node);
					}
				} else
				{
					node->nbChildren = 0;
					node->children.clear();
					FIFO.push(node);
				}	
			}
		}
	}
	return FIFO;
}



void Tree::addMutationsToNode(TreeNode* node)
{
	if (SSP->T4_MutationRate.size() == 1)
	{
		// How many mutations?
		double MRfrom = SSP->T4_MutationRate[0] * (node->left + 1);
	    double MRto   = SSP->T4_MutationRate[0] * (node->right + 1);
	    std::poisson_distribution<uint32_t> poissonDist(MRto - MRfrom);
		uint32_t nbMutations = poissonDist(GP->rngw.getRNG());

	    // Place mutations
	    for (uint32_t mutation = 0 ; mutation < nbMutations ; mutation++)
	    {
	    	// Find position
	    	std::uniform_int_distribution<uint32_t> dist(0, node->genotype.size()-1);
	        auto MutPosition = dist(GP->rngw.getRNG());

	        // toggle bit
        	node->genotype[MutPosition] = !(node->genotype[MutPosition]);
	    }
		    

	} else
	{
		// How many mutations?
		auto itFrom = SSP->T4_MutationRate.begin() + node->left;
	    auto itTo = SSP->T4_MutationRate.begin() + node->right;
	    std::poisson_distribution<uint32_t> poissonDist(*itTo - *itFrom);
		uint32_t nbMutations = poissonDist(GP->rngw.getRNG());

		// Place mutations
	    for (uint32_t mutation = 0 ; mutation < nbMutations ; mutation++)
	    {
	    	// Find position - rnd
	    	std::uniform_real_distribution<double> dist(*itFrom, *itTo);
            double rnd = dist(GP->rngw.getRNG());
            
            // Find position - binary search
            auto MutPosition = distance(itFrom,
				std::upper_bound(itFrom, itTo, rnd)
			);

	        // toggle bit
        	node->genotype[MutPosition] = !(node->genotype[MutPosition]);
	    }
	}
}

std::vector<TreeNode*> Tree::placeMutationsOnTree(std::queue<TreeNode*> FIFO, TreeNode* lastOfGenerationNode)
{
#ifdef DEBUG
	std::cout << "Enters Tree::placeMutationsOnTree\n";
#endif	
	// Explore the tree and place mutations
	std::vector<TreeNode*> nodesToDestroy;
	while (FIFO.size() > 0)
	{
		//std::cout << "FIFO.size() = " << FIFO.size() << "\n";
		TreeNode* node = FIFO.front(); FIFO.pop();
		assert(node->genotype.size() == node->right - node->left);
		assert(node->nbChildren == node->children.size());
		assert(node->nbChildren > 0);

		for (uint32_t childIndex = 0 ; childIndex < node->nbChildren ; childIndex++)
		{
			//std::cout << "childIndex = " << childIndex << "\n";
			auto& child = node->children[childIndex];
			
			assert(child->left >= node->left);
			assert(child->right <= node->right);
			// set child's genotype and free memory from node's genotype if last child
			if (childIndex == node->children.size() - 1)
			{
				// If it is last child, then I have to make sure to free the memory from the parent
				if (child->left == node->left)
				{
					// just swap. 1) it frees parents memory 2) it is faster to swap
					assert(child->genotype.size() == 0);
					child->genotype.swap(node->genotype);

					// correct right side if needed
					if (child->right < node->right)
					{
						child->genotype.resize(child->right - child->left);
					}
				} else
				{
					/*
					std::cout << "node->genotype.size() = " << node->genotype.size() << "\n";
					std::cout << "child->left = " << child->left << " child->right = " << child->right << "\n";
					*/
					
					std::vector<bool>::iterator itFrom = node->genotype.begin() + child->left - node->left;
					std::vector<bool>::iterator itTo = node->genotype.end() - node->right + child->right;
					child->genotype = std::vector<bool>(itFrom, itTo);
				}
				node->genotype.clear();
				node->genotype.shrink_to_fit();
			} else
			{
				// copy parent genotype
				/*
				std::cout << "node->genotype.size() = " << node->genotype.size() << "\n";
				std::cout << "child->left = " << child->left << " child->right = " << child->right << "\n";
				*/
				std::vector<bool>::iterator itFrom = node->genotype.begin() + child->left - node->left;
				std::vector<bool>::iterator itTo = node->genotype.end() - node->right + child->right;
				child->genotype = std::vector<bool>(itFrom, itTo);
			}
			
			// Check mutations to add
			addMutationsToNode(child);
			
			// incremenet generation if needed
			if (node == lastOfGenerationNode && childIndex == node->nbChildren - 1)
			{
				lastOfGenerationNode = child;
			}
			// add child into FIFO if not current generation
			if (child->nbChildren > 0)
				FIFO.push(child);
		}
		node->announceDeathToChildren();
		nodesToDestroy.push_back(node);
	}

#ifdef DEBUG
	std::cout << "Exits Tree::placeMutationsOnTree\n";
#endif		
	return nodesToDestroy;
}



void Tree::computeCurrentStates()
{
#ifdef DEBUG
  std::cout << "Enters in Tree::computeCurrentStates\n";
#endif	
  	//std::cout << "   Computing current states from the coalescent tree for T4 loci ... ";
	assert(!isCurrentlyBuildingHaplotype);
	assert(!currentStatesGotJustComputed);

	// initialize FIFO
	auto roots = getRootsAndSetChildren();
	std::queue<TreeNode*> FIFO;
	for (auto& root : roots)
	{	
		assert(root->nbChildren == root->children.size());
		FIFO.push(root);
	}
	TreeNode* lastOfGenerationNode = roots.back();

	// print tree
	if (file != nullptr)
	{
		printToFile(FIFO);
	}
	
	// Place mutations on tree and compute nodes to destroy
	std::vector<TreeNode*> nodesToDestroy = placeMutationsOnTree(FIFO, lastOfGenerationNode);

	if (!outputWriter.isFile(T4CoalescenceFst))
	{
		// Get rid of the tree prior to the current generation
		std::reverse(nodesToDestroy.begin(), nodesToDestroy.end());
		destroyListOfNodes(nodesToDestroy); // reverse to start with the children
		lastGenerationHaplotypes.resize(0);

		// mergeCurrentHaplotypes
		mergeCurrentHaplotypes(); // also remove parent from nodes.		
	}
		
	// know that current state has just been computed
	currentStatesGotJustComputed = true;
	//std::cout << "Done\n";
#ifdef DEBUG
  std::cout << "Exits Tree::computeCurrentStates\n";
#endif		
}



std::vector<TreeNode*> Tree::getRootsAndSetChildren() // I am assuming I can use a reference and modify the original object
{	
	TreeNode* lastOfGenerationNode = nullptr;
	std::queue<TreeNode*> FIFO = initiateFIFO(currentGenerationHaplotypes, false); // all nodes in currentGenerationHaplotypes
	std::vector<TreeNode*> roots;
	while (FIFO.size() != 0)
	{
		//std::cout << "FIFO.size() = "<< FIFO.size() << "\n";
		TreeNode* node = FIFO.front(); FIFO.pop();
		lastOfGenerationNode = node;

		if (node->parent == nullptr)
		{
			assert(node->nbChildren > 0);
			roots.push_back(node);
		} else
		{
			if (node->parent->children.size() == 0)
			{
				// Then the parent has never been added to FIFO
				FIFO.push(node->parent); 
			}
			node->parent->children.push_back(node);			
		}
	}

	std::sort( roots.begin(), roots.end() );
	roots.erase( std::unique( roots.begin(), roots.end() ), roots.end() );

	assert(roots.size() > 0 && roots.size() < SSP->TotalpatchCapacity * 2 );

	return roots;

}
	



void Tree::mergeCurrentHaplotypes()
{
	std::vector<std::vector<HaplotypeDescription>> mergedCurrentGenerationHaplotypes(currentGenerationHaplotypes.size());
	assert(currentGenerationHaplotypes.size() == GP->PatchNumber);
	for (uint32_t patch_index = 0 ; patch_index < currentGenerationHaplotypes.size() ; patch_index++)
	{
		mergedCurrentGenerationHaplotypes[patch_index].resize(currentGenerationHaplotypes[patch_index].size());
		assert(currentGenerationHaplotypes[patch_index].size() == 2 * SSP->patchSize[patch_index]);
		for (uint32_t indHaplo_index = 0 ; indHaplo_index < currentGenerationHaplotypes[patch_index].size() ; indHaplo_index++)
		{
			mergedCurrentGenerationHaplotypes[patch_index][indHaplo_index] = mergeHaplotype(currentGenerationHaplotypes[patch_index][indHaplo_index]);
		}
	}

	currentGenerationHaplotypes = mergedCurrentGenerationHaplotypes;
}

HaplotypeDescription Tree::mergeHaplotype(HaplotypeDescription& original)
{
	uint32_t lastRight = 0;
	TreeNode* newNode = new TreeNode(0, SSP->Gmap.T4_nbLoci);
	nbNodesInTree++;
	for (auto& node : original)
	{
		// security
		assert(lastRight == node->left);
		lastRight = node->right;
	
		// mergeGenotype
		newNode->genotype.insert(newNode->genotype.end(), node->genotype.begin(), node->genotype.end());
		
		// delete node
		destroyNode(node);
	}

	return {newNode};
}

std::pair<int,int> Tree::findNodeIDin(TreeNode*& node, std::vector<std::vector<HaplotypeDescription>>& listOfNodes)
{
#ifdef DEBUG
  std::cout << "Enters in Tree::findNodeIDin\n";
#endif		
	std::pair<int, int> IDnode = {-1,-1};
	assert(listOfNodes.size() == GP->PatchNumber);
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
	{
		for (uint32_t indHaplo_index = 0 ; indHaplo_index < listOfNodes[patch_index].size() ; indHaplo_index++)
		{
			for (auto& chunk : listOfNodes[patch_index][indHaplo_index])
				if (chunk == node)
				{
					// node ID found
					IDnode = {patch_index, indHaplo_index};
					goto endOfLoopCheckGoodID;
				}
		}
	}

	endOfLoopCheckGoodID:	
	if (IDnode.first < 0 || IDnode.second < 0)
	{
		std::cout << "Internal error in 'Tree::findNodeIDin'. It appears that a node can not be found where it was expected to be found\n";
		abort();
	}

	return IDnode;
}

void Tree::swapLastAndcurrentGenerationHaplotypes()
{
	//std::cout << "lastGenerationHaplotypes.size() = " << lastGenerationHaplotypes.size() << "\n";
	//std::cout << "currentGenerationHaplotypes.size() = " << lastGenerationHaplotypes.size() << "\n";
	lastGenerationHaplotypes.swap(currentGenerationHaplotypes);
	//std::cout << "lastGenerationHaplotypes.size() = " << lastGenerationHaplotypes.size() << "\n";
	//std::cout << "currentGenerationHaplotypes.size() = " << lastGenerationHaplotypes.size() << "\n";
}

void Tree::newGeneration()
{
#ifdef DEBUG
  std::cout << "Enters in Tree::newGeneration\n";
#endif
  /*
  std::cout << "TotalNbNodesCurrentGeneration = " << TotalNbNodesCurrentGeneration << "\n";
  std::cout << "TotalNbIndividualsCurrentGeneration = " << TotalNbIndividualsCurrentGeneration << "\n";
  std::cout << "SSP->T4_maxAverageNbNodesPerHaplotypeBeforeRecalculation = " << SSP->T4_maxAverageNbNodesPerHaplotypeBeforeRecalculation << "\n";
  */

	this->swapLastAndcurrentGenerationHaplotypes();
	// reset the new currentGenerationHaplotypes
	currentGenerationHaplotypes.resize(GP->PatchNumber);
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
	{
		currentGenerationHaplotypes[patch_index].resize(0);
	}
	// make sure to indicate ancestral state is not same as current state
	currentStatesGotJustComputed = false;
	TotalNbNodesCurrentGeneration = 0;
	TotalNbIndividualsCurrentGeneration = 0;
#ifdef DEBUG
  std::cout << "Exits Tree::newGeneration\n";
#endif	
}

void Tree::addChildHaplotype_setParentInfo(uint32_t parent_patch, uint32_t parent_ind_index)
{
	parent_patch_index_newHaplotypeDescription = parent_patch;
	parent_ind_index_newHaplotypeDescription   = parent_ind_index;
	assert(!isCurrentlyBuildingHaplotype);
	isCurrentlyBuildingHaplotype               = true;
	//std::cout << "parentInfo set\n";
}

void Tree::addChildHaplotype_addNode(uint32_t haploIndex, uint32_t from, uint32_t to)
{
// from included (just like left)
// to excluded   (just like right)
#ifdef DEBUG
  //std::cout << "Enters in Tree::addChildHaplotype_addNode: from = "<<from<<" to = "<<to<<"\n";
  //std::cout << "lastGenerationHaplotypes.size() = " << lastGenerationHaplotypes.size() << "\n";
  //std::cout << "parent_patch_index_newHaplotypeDescription = " << parent_patch_index_newHaplotypeDescription << "\n";
  //std::cout << "lastGenerationHaplotypes["<<parent_patch_index_newHaplotypeDescription<<"].size() = " << lastGenerationHaplotypes[parent_patch_index_newHaplotypeDescription].size() << "\n";
  //std::cout << "parent_ind_index_newHaplotypeDescription = " << parent_ind_index_newHaplotypeDescription << "\n";
  //std::cout << "haploIndex = " << haploIndex << "\n";
  assert(lastGenerationHaplotypes.size() > parent_patch_index_newHaplotypeDescription);
  assert(lastGenerationHaplotypes[parent_patch_index_newHaplotypeDescription].size() > parent_ind_index_newHaplotypeDescription * 2 + haploIndex);
#endif


	HaplotypeDescription& parentHaplotype = lastGenerationHaplotypes[parent_patch_index_newHaplotypeDescription][parent_ind_index_newHaplotypeDescription * 2 + haploIndex];

	for (uint32_t ParentHaplotypeChunkIndex = 0;ParentHaplotypeChunkIndex < parentHaplotype.size() ; ParentHaplotypeChunkIndex++)
	{
		TreeNode* chunk = parentHaplotype[ParentHaplotypeChunkIndex];
		
		/*
		std::cout << "ParentHaplotypeChunkIndex = " << ParentHaplotypeChunkIndex << " --- ";
		std::cout << "chunk->left = " << chunk->left << "  chunk->right = "<< chunk->right<< " --- ";
		std::cout << "from = " << from << "  to = "<< to << "\n";
		assert(to <= SSP->Gmap.T4_nbLoci);
		*/
		

		assert(chunk->left < to); // Otherwise, it should have been in a previous chunk and has been missed somehow

		
		if (chunk->right <= from)
		{
			// Then it is not in chunk
			continue;
		} else
		{
			// Then it is in this chunk

			TreeNode* newNode;
			if (chunk->right >= to)
			{
				// Then it ends in this chunk
				if (chunk->left <= from)
				{
					// Then it started in this chunk
					newNode = new TreeNode(chunk, from, to);	
				} else
				{
					// Then it started in a previous chunk
					newNode = new TreeNode(chunk, chunk->left, to);
				}
				this->addNode(newNode); // can't move this addNode because opf the break; statement that follows
				break;
			} else
			{
				// It ends after this chunk
				if (chunk->left <= from)
				{
					// Then it started in this chunk
					newNode = new TreeNode(chunk, from, chunk->right);
				} else
				{
					// Then it started in a previous chunk
					newNode = new TreeNode(chunk, chunk->left, chunk->right);
				}
				this->addNode(newNode);
				continue;
			}
		}
	}
	
	/*
	for (TreeNode*& chunk : parentHaplotype)
	{
		uint32_t newleft = chunk->left;
		uint32_t newright;
		for (auto& breakpoint : T4breakpoints)
		{
			if (chunk->left <= breakpoint && breakpoint < chunk->right )
			{
				newright = breakpoint;
				TreeNode* newNode = new TreeNode(chunk, newleft, newright);
				this->addNode(newNode);
				newleft = breakpoint + 1;
			}
		}
		newright = chunk->right;
		if (newleft != newright + 1)
		{
			TreeNode* newNode = new TreeNode(chunk, newleft, newright);
			this->addNode(newNode);
		}
	}
	this->addChildHaplotype_finished(child_patch_index);
	*/
#ifdef DEBUG
  std::cout << "Exits Tree::addChildHaplotype_addNode\n";
#endif			
}



void Tree::pruneDeadLineages()
{
	if (!currentStatesGotJustComputed)
	{
		auto FIFO = initiateFIFO(lastGenerationHaplotypes,true);
		pruneFromFIFO(FIFO);

		if (
			!outputWriter.isFile(T4CoalescenceFst)
			&& 
			(TotalNbIndividualsCurrentGeneration > 0 && (double) TotalNbNodesCurrentGeneration / (double) TotalNbIndividualsCurrentGeneration > SSP->T4_maxAverageNbNodesPerHaplotypeBeforeRecalculation)
		)
		{
			computeCurrentStates();
		}
	}
}



std::vector<std::vector<double>> Tree::getCurrentStates_frequencies()
{
	// compute if needed
	if (!currentStatesGotJustComputed)
	{
		computeCurrentStates();	
	}

	std::vector<std::vector<double>> r(GP->PatchNumber);
	for (uint32_t patch_index =0 ; patch_index < GP->PatchNumber ; patch_index++)
	{
		// Sum up the trues
		r[patch_index].resize(SSP->Gmap.T4_nbLoci);
		for (uint32_t ind_index =0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)	
		{
			for (uint32_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
			{
				assert(currentGenerationHaplotypes[patch_index][ind_index * 2 + haplo_index].size() == 1);
				assert(currentGenerationHaplotypes[patch_index][ind_index * 2 + haplo_index][0]->genotype.size() == SSP->Gmap.T4_nbLoci);

				for (uint32_t locus = 0 ; locus < SSP->Gmap.T4_nbLoci ; locus++)
				{
					if (currentGenerationHaplotypes[patch_index][ind_index * 2 + haplo_index][0]->genotype[locus])
						r[patch_index][locus]++;
				}
			}
		}

		// divide
		for (uint32_t locus = 0 ; locus < SSP->Gmap.T4_nbLoci ; locus++)       
	  	{
	    	r[patch_index][locus] /= 2 * SSP->patchSize[patch_index];
	    	assert(r[patch_index][locus] >= 0.0 && r[patch_index][locus] <= 1.0);
	  	}
	}

	
	return r;	
}


std::vector<std::vector<std::vector<std::vector<bool>>>> Tree::getCurrentStates()
{
	// compute if needed
	if (!currentStatesGotJustComputed)
	{
		computeCurrentStates();	
	}

	// Go from nodes to a nice looking vector
	std::vector<std::vector<std::vector<std::vector<bool>>>> r(GP->PatchNumber);
	for (uint32_t patch_index =0 ; patch_index < GP->PatchNumber ; patch_index++)
	{
		r[patch_index].resize(SSP->patchSize[patch_index]);
		for (uint32_t ind_index =0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)	
		{
			r[patch_index][ind_index].resize(2);
			for (uint32_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
			{
				assert(currentGenerationHaplotypes[patch_index][ind_index * 2 + haplo_index].size() == 1);
				assert(currentGenerationHaplotypes[patch_index][ind_index * 2 + haplo_index][0]->genotype.size() == SSP->Gmap.T4_nbLoci);
				r[patch_index][ind_index][haplo_index] = currentGenerationHaplotypes[patch_index][ind_index * 2 + haplo_index][0]->genotype;
			}
		}
	}
	
	return r;
}


