


void T4TreeRec::initialize(Pop& pop)
{
	/////////////////
	// Reserve RAM //
	/////////////////

	{
		// How much RAM to reserve?
		unsigned long long maxEverNbNodesProduced = 0;
		assert(GP->__GenerationChange.size() == SSP->__patchCapacity.size());
		for (size_t generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; ++generation_index )
		{
			int nbGens;
			if (generation_index == GP->__GenerationChange.size()-1)
			{
				nbGens = GP->nbGenerations - GP->__GenerationChange[generation_index];
			} else
			{
				nbGens = GP->__GenerationChange[generation_index+1] - GP->__GenerationChange[generation_index];
			}
			assert(nbGens >= 0);


			int totalPatchCapacity = 0;
			for (auto& pc : SSP->__patchCapacity[generation_index])
			{
				totalPatchCapacity += pc;
			}
			
			maxEverNbNodesProduced += (unsigned long long) totalPatchCapacity * (unsigned long long) nbGens;
		}

		auto nbNodesToReserveFor = std::min(maxEverNbNodesProduced, (unsigned long long) SSP->T4_maxNbEdges);

		// Reserve RAM
		// nodes.data.reserve(nbNodesToReserveFor); it is a deque. No need to reserve
		edges.data.reserve(nbNodesToReserveFor * (1 + SSP->TotalRecombinationRate));
	}
	

/*	/////////////////////////////////
	// Add Ancestor and set pop ID //
	/////////////////////////////////

	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber; ++patch_index )
	{
		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index )
		{
			for (size_t haplo_index = 0 ; haplo_index < 2; ++haplo_index )
			{
				NodeGenetics ancestorGenetics; // initialize at 0 for all loci
				auto T4ID = addAncestor(ancestorGenetics); // Just one node is enough as they are all clones

				pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).T4ID = T4ID;
			}
		}
	}
*/


	//////////////////
	// Add Ancestor //
	//////////////////

	NodeGenetics ancestorGenetics; // initialize at 0 for all loci
	auto T4ID = addAncestor(ancestorGenetics); // Just one node is enough as they are all clones

	/////////////////////////
	// Set population T4ID //
	/////////////////////////

	setPopToUniqueID(pop, T4ID);
}

void T4TreeRec::setPopToUniqueID(Pop& pop, size_t T4ID)
{
	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber; ++patch_index )
	{
		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index )	
		{
			auto& haplo0 = pop.getPatch(patch_index).getInd(ind_index).getHaplo(0);
			auto& haplo1 = pop.getPatch(patch_index).getInd(ind_index).getHaplo(1);

			haplo0.T4ID = T4ID;
			haplo1.T4ID = T4ID;
		}
	}
}


void T4TreeRec::printInfoForDebug(Pop& pop) 
{
	// Nb nodes per generation
	std::vector<size_t> nbNodesPerGeneration(GP->CurrentGeneration+1, 0);
	for (size_t u = 0 ; u < nodes.size() ; ++u )
	{
		assert(nodes[u].birth >= 0 && nodes[u].birth < nbNodesPerGeneration.size());
		++nbNodesPerGeneration[nodes[u].birth];
	}

	size_t sum = 0;
	for (size_t gen = 0 ; gen < nbNodesPerGeneration.size() ; ++gen)
	{
		std::cout << "Nb nodes in generation: " <<gen << " is " << nbNodesPerGeneration[gen] << "\n";
		sum += nbNodesPerGeneration[gen];
	}
	assert(sum == nodes.size());


	// What fraction of haplotypes have parented at each generation
	std::unordered_map<size_t, std::vector<size_t>> hash;
	for (size_t edge_index = 0 ; edge_index < edges.size() ; ++edge_index)
	{
		hash[edges[edge_index].parent].push_back(edge_index);
	}

	std::vector<size_t> nbHaplosThatReproducedAtEachGeneration(GP->CurrentGeneration+1, 0);
	for (size_t u = 0 ; u < nodes.size() ; ++u )
	{
		if (hash[u].size())
		{
			++nbHaplosThatReproducedAtEachGeneration[nodes[u].birth];
		}
	}

	assert(nbHaplosThatReproducedAtEachGeneration.size() == nbNodesPerGeneration.size());
	for (size_t gen = 0 ; gen < nbNodesPerGeneration.size() ; ++gen)
	{
		assert(nbHaplosThatReproducedAtEachGeneration.size() > gen);
		assert(nbNodesPerGeneration.size() > gen);

		std::cout << "Fractions of haplodes that reproduced at generation: " <<gen << " is " << (double)nbHaplosThatReproducedAtEachGeneration[gen] / (double)nbNodesPerGeneration[gen] << "\n";
	}

	

	// Pop IDs
	std::vector<size_t> IDtable;
	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			for (size_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
			{
				auto ID = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).T4ID;
				assert(ID >= 0 && ID <nodes.size());
				if (ID >= IDtable.size())
				{
					IDtable.resize(ID+1, 0);
				}

				++IDtable[ID];
			}
		}
	}

	/*for (size_t ID = 0 ; ID < IDtable.size() ; ++ID)
	{
		std::cout << "popID: " << ID << ": " << IDtable[ID] << "\n";
	}	*/

}


/*class mapT4TreeStructure
{
	std::vector<size_t> birthTimeBoundaries;

	void compute(NodeTable& nodes)
	{
		birthTimeBoundaries.resize(0);

		for (size_t nodeID = 0 ; nodeID < nodes.size() ; ++nodeID)
		{

		}
	}

	size_t beginGeneration()
	{

	}

	size_t endGeneration()
	{

	}


};
*/


std::vector<std::vector<std::vector<size_t>>> T4TreeRec::placeMutations(Pop& pop, bool isNeedSimplify) // Pop is used to access haplotypes T4ID to do a good matching and to reset them
{
/*
	It computes the current states and set all current individuals as ancestors (unless if there are clones, then they share the same T4ID)
*/

	//////////////////////////////////////////////////////////////////////////////
	//// If I don't need to place mutations but just return the current state ////
	//////////////////////////////////////////////////////////////////////////////


	if (nodes.size() == ancestorsGenetics.size())
	{
		assert(ancestorsGenetics.size() == ancestorsID.size());

		std::vector<std::vector<std::vector<size_t>>> ret(GP->PatchNumber);
		for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
		{
			auto& patch = pop.getPatch(patch_index);
			ret[patch_index].reserve(2*SSP->patchSize[patch_index]);

			assert(patch.getpatchCapacity() >= SSP->patchSize[patch_index]);
			for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
			{
				auto& ind = patch.getInd(ind_index);

				// Get haplotypes
				auto& haplo0 = ind.getHaplo(0);
				auto& haplo1 = ind.getHaplo(1);

				// Assert
				assert(ancestorsGenetics.size() > haplo0.T4ID);
				assert(ancestorsGenetics.size() > haplo1.T4ID);

				// Build ret (object to return)
				ret[patch_index].push_back(ancestorsGenetics[haplo0.T4ID].mutations);
				ret[patch_index].push_back(ancestorsGenetics[haplo1.T4ID].mutations);
			}	
		}


		return ret;
	}



	assert(SSP->T4_MutationRate.size() == 1 || SSP->T4_MutationRate.size() == SSP->T4_nbLoci);

	////////////////////////
	//// Simplify first ////
	////////////////////////

	//printInfoForDebug(pop);

	if (isNeedSimplify && edges.size()) simplify(pop);
	
	//printInfoForDebug(pop);
	

	/////////////////////////////////
	//// Hash table to offspring ////
	/////////////////////////////////

	// Note I can't reuse the hash table from simplify because simplify simplified the tree since simplify's hash table has been created.

	std::unordered_map<size_t, std::vector<size_t>> hash;
	for (size_t edge_index = 0 ; edge_index < edges.size() ; ++edge_index)
	{
		assert(edges[edge_index].parent >= 0 && edges[edge_index].parent < nodes.size());
		hash[edges[edge_index].parent].push_back(edge_index);
	}

	//////////////////////////////////
	//// initialize Tree genetics ////
	//////////////////////////////////

	std::vector<NodeGenetics> treeGenetics(nodes.size());


	//////////////////////////////////////
	//// Take the ancestral variation ////
	//////////////////////////////////////

	assert(nodes.size());
	assert(ancestorsGenetics.size() == ancestorsID.size());
	for (size_t ancestor_index = 0 ; ancestor_index < ancestorsID.size() ; ++ancestor_index)
	{
		assert(treeGenetics.size() > ancestorsID[ancestor_index]);
		treeGenetics[ancestorsID[ancestor_index]] = ancestorsGenetics[ancestor_index];
		//treeGenetics[ancestorsID[ancestor_index]].whatHasBeenSetYet.push({0, (size_t) SSP->T4_nbLoci,-1});
		//std::cout << "Setting ID " << ancestorsID[ancestor_index] << " in treeGenetics. This individual has birth = "<< ancestorsGenetics[ancestor_index].birth <<"\n";
	}


	////////////////////////////////////////////////////
	//// Set entire treeGenetics and set new 'this' ////
	////////////////////////////////////////////////////

	//std::cout << "About to set entire table\n";

	NodeTable oldNodes;
	oldNodes.swap(nodes); // nodes is now empty
	ancestorsGenetics.resize(0);
	ancestorsID.resize(0);
	isAncestor.resize(0);
	//std::cout << "oldNodes.size() = " << oldNodes.size() << "\n";

	size_t nbIndsBornInCurrentGeneration = 0; // Just a security
	std::vector<size_t> oldToNewIDMap(oldNodes.size());
	for (size_t u = 0 ; u < oldNodes.size() ; ++u )
	{
		// Loop through children
		auto& edges_u_isParent = hash[u];
		for (auto& edge_index : edges_u_isParent)
		{
			auto& edge = edges[edge_index];
			assert(edge.parent == u);
			//std::cout << "Using IDs " << edge.parent << " as parent and ID "<<edge.child<<" as child\n";
			assert(edge.parent < oldNodes.size());
			assert(edge.child < oldNodes.size());
			//std::cout << "About to set birth .. ";
			//std::cout << ".. to " << oldNodes[edge.child].birth << "\n";
			treeGenetics[edge.child].birth = oldNodes[edge.child].birth;
			//std::cout << "u = "<<u<<". This individual has birth = "<<treeGenetics[edge.parent].birth<<"\n";

			treeGenetics[edge.child].propagateMutationsForSegment(treeGenetics[edge.parent], edge.left, edge.right);
			//std::cout << "Propagate mutations from " << edge.parent << " to " << edge.child << "\n";
		}

		
		// Create new tree and new ancestors
		if (oldNodes[u].birth == GP->CurrentGeneration)
		{

			assert(treeGenetics[u].birth == oldNodes[u].birth);
			++nbIndsBornInCurrentGeneration; // just a security

			assert(edges_u_isParent.size() == 0);

			int cloneAncestorID = searchForClonesAmongAncestors(treeGenetics[u]);
			size_t v;
			if (cloneAncestorID == -1) // no clone found
			{
				v = addAncestor(treeGenetics[u]); // Ancestor added into nodes
				assert(ancestorsGenetics.size() > v);
				//std::cout << "Asserted ancestor ID " << v << "\n";
				assert(v == nodes.size()-1);
				/*
				std::cout << "Added ancestor oldID " << u << ", newID " << v << " with mutations:";
				for (auto& mut : treeGenetics[u].mutations)
					std::cout << mut << " ";
				std::cout << "\n";
				*/
			} else
			{
				v = cloneAncestorID;
				/*
				std::cout << "Ancestor oldID " << u << " is a clone with newID " << v << " with mutations:";
				for (auto& mut : treeGenetics[u].mutations)
					std::cout << mut << " ";
				std::cout << "\n";
				*/
			}

			oldToNewIDMap[u] = v;
		}

		// Remove RAM that won't be used anymore
		treeGenetics[u].clear();
		treeGenetics[u].shrink_to_fit();
	}
	
	assert(nbIndsBornInCurrentGeneration == SSP->TotalpatchSize*2);
	edges.clear(); // No edges in the new tree
	assert(nodes.size());
	assert(nodes.size() <= 2*SSP->TotalpatchCapacity);
	assert(ancestorsID.size() == nodes.size());


	////////////////////////////////////////////////////////
	//// Gather data in convenient format and reset IDs ////
	////////////////////////////////////////////////////////

	//std::cout << "About to gather data\n";

	std::vector<std::vector<std::vector<size_t>>> ret(GP->PatchNumber);
	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);
		ret[patch_index].reserve(2*SSP->patchSize[patch_index]);

		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			auto& ind = patch.getInd(ind_index);

			// Get haplotypes
			auto& haplo0 = ind.getHaplo(0);
			auto& haplo1 = ind.getHaplo(1);

			// Assertion
			assert(oldToNewIDMap[haplo0.T4ID] < 2*SSP->patchSize[patch_index]);
			assert(oldToNewIDMap[haplo1.T4ID] < 2*SSP->patchSize[patch_index]);

			// reset the IDs
			haplo0.T4ID = oldToNewIDMap[haplo0.T4ID];
			haplo1.T4ID = oldToNewIDMap[haplo1.T4ID];

			// Assert
			assert(ancestorsGenetics.size() > haplo0.T4ID);
			assert(ancestorsGenetics.size() > haplo1.T4ID);

			// Build ret (object to return)
			//std::cout << "treeGenetics["<<haplo0.T4ID<<"].mutations: "; printVector(treeGenetics[haplo0.T4ID].mutations);
			//std::cout << "treeGenetics["<<haplo1.T4ID<<"].mutations: "; printVector(treeGenetics[haplo1.T4ID].mutations);
			assert(ancestorsID[haplo0.T4ID] == haplo0.T4ID);
			assert(ancestorsID[haplo1.T4ID] == haplo1.T4ID);
			//std::cout << "Asserting ancestor ID " << haplo0.T4ID << "...";
			//std::cout << "Done\n";
			//std::cout << "Asserting ancestor ID " << haplo1.T4ID << "...";
			//std::cout << "Done\n";
			ret[patch_index].push_back(ancestorsGenetics[haplo0.T4ID].mutations);
			ret[patch_index].push_back(ancestorsGenetics[haplo1.T4ID].mutations);
		}	
	}
	//std::cout << "Finish gathering data\n";


	return ret;
}


int T4TreeRec::searchForClonesAmongAncestors(const NodeGenetics& input) const
{
	assert(ancestorsGenetics.size() == ancestorsID.size());

	for (size_t ancestor_index = 0 ; ancestor_index < ancestorsGenetics.size() ; ++ancestor_index)
	{
		const auto& ancestor = ancestorsGenetics[ancestor_index];
		
		// test if ancestor is clone of input
		if (input.mutations.size() != ancestor.mutations.size())
		{
			continue;
		} else
		{
			bool areClones = true;
			for (size_t mut_index = 0 ; mut_index < input.mutations.size() ; ++mut_index )
			{
				if (input.mutations[mut_index] != ancestor.mutations[mut_index])
				{
					areClones = false;
					break;
				}
			}
			if (!areClones) continue;
		}

		// If it gets there it is because it did not goto nextAncestor and therefore they are clone

		/*
		std::cout << "Clone found. Input muts: ";
		for (auto& mut : input.mutations) std::cout << mut << " ";
		std::cout << " Ancestor muts: ";
		for (auto& mut : ancestor.mutations) std::cout << mut << " ";
		std::cout << "\n";
		*/

		return ancestorsID[ancestor_index];
	}
	//std::cout << "Clone not found\n";
	return -1;
}

template<typename INT>
void NodeGenetics::mutateLocus(INT MutPosition)
{
	auto position = std::lower_bound(mutations.begin(), mutations.end(), MutPosition);

    if (position == mutations.end())
    {
        // not found
        mutations.push_back(MutPosition);  
    } else
    {
        if ( MutPosition == (*position))
        {
            // found
            mutations.erase(position);
        } else
        {
            // not found
            mutations.insert(position, MutPosition);
        }
    }
}

/*
void NodeGenetics::assertWasFullySet() const
{
	auto Q = whatHasBeenSetYet; // copy

	assert(Q.size());
	auto firstSeg = Q.top(); Q.pop();
	assert(firstSeg.left == 0);
	auto nextLeft = firstSeg.right;

	while (Q.size())
	{
		auto& S = Q.top(); Q.pop();
		assert(nextLeft == S.left);
		nextLeft = S.right;
	}
	assert(nextLeft == SSP->T4_nbLoci);
}
*/

void NodeGenetics::propagateMutationsForSegment(const NodeGenetics& parent, const size_t left, const size_t right)
{
	assert(left < right);
	assert(this->birth > parent.birth);
	auto nbGenerationsInBetween = this->birth - parent.birth;

	//whatHasBeenSetYet.push({left, right, -1});


	////////////////////////////////
	// Set segment as from parent //
	////////////////////////////////

	if (parent.mutations.size())
	{
		// Get iterators
		auto insertFrom = std::lower_bound(parent.mutations.cbegin(), parent.mutations.cend(), left);    // lower only works
		auto insertTo = std::lower_bound(insertFrom, parent.mutations.cend(), right);                   // lower only works
		auto startOfInsertion = std::lower_bound(this->mutations.begin(), this->mutations.end(), left); // lower or upper. Both should work

		// security -> Only one parent node can affect the a given segment
			
		if (startOfInsertion != this->mutations.end()) 
		{
			
			// if something comes after the point of insertion then the next value must be greater or equal to right
			assert( *startOfInsertion >= left );

			if (startOfInsertion != (this->mutations.end() - 1))
			{
				assert( *(startOfInsertion+1) >= right );
			}
		}

		// insert
		if (insertFrom != insertTo)
			this->mutations.insert(startOfInsertion, insertFrom, insertTo);
/*
		if (mutations.size()>1)
		{
			size_t prev = std::numeric_limits<size_t>::max();
			for (auto& m : mutations)
			{
				if (prev != std::numeric_limits<size_t>::max() && prev >= m)
				{
					std::cout << "----> JUST BEFORE ABORTING\nFrom parent " << pID << " to child " << cID << "\n";
					std::cout << "left = " << left << "\n";
					std::cout << "right = " << right << "\n";
					std::cout << "child mutations:";
					for (auto& m : mutations) std::cout << m << " ";
					std::cout << "\n";
					std::cout << "parent mutations:";
					for (auto& m : parent.mutations) std::cout << m << " ";
					std::cout << "\nstartOfInsertion - mutations.begin() = " << startOfInsertion - mutations.begin() << "\n";
					std::cout << "insertFrom - parent.mutations.begin() = " << insertFrom - parent.mutations.begin() << "\n";
					std::cout << "insertTo - parent.mutations.begin() = " << insertTo - parent.mutations.begin() << "\n";
					std::cout << "------------\n";
					abort();
				}
				prev = m;
			}
		}*/
	}
	

	////////////
	// Mutate //
	////////////
	int nbMuts = SSP->geneticSampler.get_T4_nbMuts(nbGenerationsInBetween, left, right);
	//std::cout << nbMuts << " muts: ";
    
    for (int i = 0 ; i < nbMuts ; i++)
    {
    	// Where does the ith mutation happen?
        auto MutPosition = SSP->geneticSampler.get_T4_mutationPosition(left, right);
        /*if (!(MutPosition>= left && MutPosition < right))
        {
        	std::cout << "right = " << right << "\n";
        	std::cout << "left = " << left << "\n";
        	std::cout << "MutPosition = " << MutPosition << "\n";
        }*/
        //std::cout << MutPosition << " ";
        assert(MutPosition>= left && MutPosition < right);

        // Make the mutation
        this->mutateLocus(MutPosition);
    }
    //std::cout <<"\n";

    /*
    std::cout << "Childs genetics: ";
    for (auto& m : child.mutations) std::cout << m << " ";
    std::cout << "\n"; 
	*/
}

size_t T4TreeRec::addAncestor(NodeGenetics& newAncestor)
{
	// Add into nodes
	auto offNodeID = nodes.addNode(GP->CurrentGeneration);

	// Set ancestor genetics
	ancestorsGenetics.push_back({newAncestor.mutations, GP->CurrentGeneration});
	//ancestorsGenetics.back().whatHasBeenSetYet = newAncestor.whatHasBeenSetYet;
	ancestorsID.push_back(offNodeID);
	if (isAncestor.size() <= offNodeID)
	{
		isAncestor.resize(offNodeID+1, false);
	}
	isAncestor[offNodeID] = true;

	//std::cout << "ancestorsID.size() = " << ancestorsID.size() << "\n";

	return offNodeID;
}

size_t T4TreeRec::addHaplotype(std::vector<int>& RP, std::pair<size_t, size_t> p)
{
	/*
		Haplotypes need to know their nodeID
	*/


	if (p.first == p.second)
	{
		RP.resize(1);
		RP[0] = INT_MAX;
	}

	//std::cout << "p.first = " << p.first << "\n";
	//std::cout << "nodes.size() = " << nodes.size() << "\n";
	assert(nodes.size() > p.first);
	assert(RP.back() == INT_MAX);


	
	auto offNodeID = nodes.addNode(GP->CurrentGeneration);
	
	if (RP.size() > 1) // This if else is just for performance
	{
		size_t l = 0;
		bool isPfirst = true;
		for (size_t DNAsegment = 0 ; DNAsegment < RP.size() ; ++DNAsegment)
		{
			size_t r = RP[DNAsegment];
			if (r == INT_MAX)
			{
				r = SSP->T4_nbLoci;
			} else
			{
				r = SSP->FromLocusToTXLocus[r].T4; 
			}
			edges.addEdge(l, r, isPfirst ? p.first : p.second, offNodeID);
			l = r;
			isPfirst = !isPfirst;
		}
	} else
	{
		edges.addEdge(0, SSP->T4_nbLoci, p.first, offNodeID);
	}

	//std::cout << "offNodeID = " << offNodeID << "\n";

	return offNodeID;
}

void T4TreeRec::simplify_ifNeeded(Pop& pop)
{
	if (
		edges.size() > SSP->T4_maxNbEdges &&
		GP->CurrentGeneration - lastGenerationSimplified > SSP->T4_minimumNbGenerationsBetweeSimplification &&
		GP->CurrentGeneration - GP->startAtGeneration > SSP->T4_minimumNbGenerationsBetweeSimplification &&
		GP->nbGenerations - GP->CurrentGeneration > SSP->T4_minimumNbGenerationsBetweeSimplification
	)
	{
		//std::cout << "Before simplification: edges.size() = "<<edges.size()<<"\n";
		//this->print();
		//auto nbEdgesBefore = edges.size();
		simplify(pop);
		//std::cout << "Simplification went from "<<nbEdgesBefore<<" edges to "<< edges.size()<<" edges.\n";
	}
}


/*
T4TreeRec::simplify_transferChildSegments_addToParentSegments(std::vector<Segment>& parentSegments, size_t left, size_t right, size_t newID)
{
	if (parentSegments.size() && parentSegments.back().right == left && parentSegments.back().newID === newID)
	{
		parentSegments.back().right = right;
	} else
	{
		parentSegments.push_back(Segment(left, right, newID));
	}
}
*/


void T4TreeRec::simplify_transferChildSegments(
	std::vector<Segment>& parentSegments,
	std::vector<Segment>& childSegments,
	size_t left,
	size_t right
)
{
	/*
		Assuming sorted segments in inputs
		This function makes plenty of little segments
	*/
	std::cout << "childSegment:\n";
	printVector(childSegments);
	std::cout << "parenttSegment:\n";
	printVector(parentSegments);

	assert(left < right);

	std::vector<Segment> addToParentSegments;

	std::cout << "checking if sorted at 723\n";
	assert(std::is_sorted(childSegments.begin(), childSegments.end()));

	std::vector<Segment> newChildSegments;

	for (size_t childSegmentIndex = 0 ; childSegmentIndex < childSegments.size() ;  ++childSegmentIndex )
	{
		auto& childSegment = childSegments[childSegmentIndex];
		std::cout << "checking accessing segment values...\n";
		std::cout << childSegment.left << " " << childSegment.right << " " << childSegment.newID << "\n";
		std::cout << "Done\n";

		if (childSegment.left < right && childSegment.right > left)
		{
			if (childSegment.right > right)
			{
				if (childSegment.left < left)
				{
					/*
					Need to duplicate child segment
					Edge truncated on both sides
						Edge:	              [---]
						childSegment:       [-------]
					*/
					
					newChildSegments.push_back({childSegment.left, left, childSegment.newID});
					newChildSegments.push_back({right, childSegment.right, childSegment.newID});
					
					addToParentSegments.push_back({left, right, childSegment.newID});

				} else
				{
					/*
					truncate left side of child segment
						Edge:	          [-------]
						childSegment:       [-------]
					*/
					if (right != childSegment.right)
						newChildSegments.push_back({right, childSegment.right, childSegment.newID});
					addToParentSegments.push_back({childSegment.left, right, childSegment.newID});
				}
			} else
			{
				if (childSegment.left < left)
				{
					/*
					truncate right side of child segment
						Edge:	               [-------]
						childSegment:       [-------]
					*/

					newChildSegments.push_back({childSegment.left, left, childSegment.newID});
					addToParentSegments.push_back({left, childSegment.right, childSegment.newID});
				} else
				{
					/*
					No more child segment. Everything is given ot parent
						Edge:	          [-----------]
						childSegment:       [-------]
					*/
					assert(childSegment.left != childSegment.right);
					addToParentSegments.push_back({childSegment.left, childSegment.right, childSegment.newID});
				}
			}
		} else
		{
			// Whole child segment unchanged. Nothing is given to parent
			newChildSegments.push_back(childSegment);
		}
	}

	std::cout << "line 797\n";	
	newChildSegments.swap(childSegments);
	std::cout << "line 799\n";	
	if (parentSegments.size())
	{
		if (addToParentSegments.size())
		{
			std::cout << "line 802\n";	
			std::vector<Segment> x;
			std::merge(
				parentSegments.begin(),
				parentSegments.end(),
				addToParentSegments.begin(),
				addToParentSegments.end(),
				x.begin()
			);
			std::cout << "line 811\n";	
			parentSegments.swap(x);
			std::cout << "line 813\n";	
		}
	} else
	{
		std::cout << "line 816\n";	
		parentSegments.swap(addToParentSegments);
		std::cout << "line 816\n";	
	}

	std::cout << "checking if sorted at 790\n";
	assert(std::is_sorted(childSegments.begin(), childSegments.end()));
	assert(std::is_sorted(parentSegments.begin(), parentSegments.end()));

}

void T4TreeRec::simplify_mergeSegmentsGoingThroughNode(std::vector<Segment>& segments)
{
	if (segments.size()>1)
	{
		std::cout << "checking if sorted at 799\n";
		assert(std::is_sorted(segments.begin(), segments.end()));

		
		for (std::vector<Segment>::iterator it = segments.begin() ; it < (segments.end() - 1) ; )
		{
			Segment& first  = *it;
			Segment& second = *(it+1);
			if (first.newID == second.newID && first.right > second.left)
			{
				// intersection
				first.right = second.right;
				segments.erase(it); // should not reallocate memory
			} else
			{
				++it;
			}
		}
	}
}

void T4TreeRec::simplify(Pop& pop)
{
	lastGenerationSimplified = GP->CurrentGeneration;
	auto& Ei = edges;
	auto& Ni = nodes;
	NodeTable No;
	EdgeTable Eo;


	
	{
		///////////////////////////////////////////////////////////////////////
		// Add latest generation and set map of IDs for this last generation //
		///////////////////////////////////////////////////////////////////////

		LastGenerationIDmap IDmap (Ni.size(), SSP->TotalpatchSize);
		size_t nbIndsFoundInGeneration = 0;
		for (size_t u = Ni.size()-1 ; u >= 0 && Ni[u].birth == Ni.back().birth ; --u)
		{
			++nbIndsFoundInGeneration;
			assert(nbIndsFoundInGeneration <= SSP->TotalpatchSize*2); 
			IDmap.add(u, No.addNode(Ni[u]));
		}

		//////////////////////////
		// Reset T4flags in Pop //
		//////////////////////////

		for (size_t patch_index = 0 ; patch_index < pop.getNbPatches() ; ++patch_index)
		{
			auto& patch = pop.getPatch(patch_index);
			for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
			{
				auto& ind = patch.getInd(ind_index);
				auto& haplo0 = ind.getHaplo(0);
				auto& haplo1 = ind.getHaplo(1);

				haplo0.T4ID = IDmap.getNewID(haplo0.T4ID);
				haplo1.T4ID = IDmap.getNewID(haplo1.T4ID);
			}	
		}
	}
	


	///////////////////////////////////////////////////////////////////
	// Set up the segments of the children of the current generation //
	///////////////////////////////////////////////////////////////////

	std::vector<std::vector<Segment>> listSegments(Ni.size());
	int childGeneration = Ni.back().birth;
	assert(childGeneration == GP->CurrentGeneration);
	int edge_index = Ei.size()-1;
	for ( int nodeID = Ni.size()-1 ; nodeID >= 0 && Ni[nodeID].birth == childGeneration ; --nodeID)
	{
		listSegments[nodeID].push_back(Segment(0, SSP->T4_nbLoci, -1));
	}



	///////////////////
	// New Ancestors //
	///////////////////

	std::vector<NodeGenetics> new_ancestorsGenetics;
	std::vector<size_t> new_ancestorsID;


	//////////////
	// Big loop //
	//////////////

	int Ni_indexForParents = Ni.size()-1;
	while (edge_index >= 0)
	{
		auto parentGeneration = childGeneration-1;

		////////////////////////////////////////////////////////
		// transfer segments from children of this generation //
		////////////////////////////////////////////////////////

		// Ensures all direct parents of generation parentGeneration know all of their segments
		while ( edge_index >= 0 && Ni[Ei[edge_index].child].birth == childGeneration)
		{
			auto parentID = Ei[edge_index].parent;
			auto childID = Ei[edge_index].child;
			assert(parentID != childID);
			this->simplify_transferChildSegments(listSegments[parentID], listSegments[childID], Ei[edge_index].left, Ei[edge_index].right);
			--edge_index;
		}


		///////////////////////////////////////////////////////
		// Reset parents segments and create nodes and edges //
		///////////////////////////////////////////////////////

		// resolve segment for each parent
		for ( ; Ni_indexForParents >= 0 && Ni[Ni_indexForParents].birth==parentGeneration ; --Ni_indexForParents)
		{
			/* INFO: A segment is an information about what a child needs to get from its parents */

			// just another name for convenience
			auto parentID = Ni_indexForParents; 
		
			
			auto& segmentsForChildren = listSegments[parentID];
			std::vector<Segment> copy_segmentsForChildren;
			if (isAncestor.size() > parentID && isAncestor[parentID]) copy_segmentsForChildren = segmentsForChildren;
			std::vector<Segment> segmentsForParents; // Those are the segments will be given to the parents of parent
			int newParentID = -1;

			while (segmentsForChildren.size())
			{
				auto segmentForChild = segmentsForChildren.back(); segmentsForChildren.pop_back();
				auto l = segmentForChild.left;
				auto r = segmentForChild.right;

				// While I find intersecting segments
				bool doesSpecificSegmentNeedToGoThroughNodeOfCurrentParent = false;
				while (segmentsForChildren.size() && segmentForChild.right > segmentsForChildren.back().left)
				{
					// There is an intersection -> A node will be added if not done yet
					if (newParentID == -1) // if not added yet
					{
						newParentID = No.addNode(Ni[parentID]);
					}
					doesSpecificSegmentNeedToGoThroughNodeOfCurrentParent = true;

					// Record edge
					Eo.addEdge(segmentForChild.left, segmentForChild.right, newParentID, segmentForChild.newID);

					// get to next segment
					r = segmentForChild.right;
					segmentForChild = segmentsForChildren.back();
					segmentsForChildren.pop_back();

				}

				// Create a segment for parents of parent
				if (doesSpecificSegmentNeedToGoThroughNodeOfCurrentParent)
				{
					segmentsForParents.push_back(Segment(l, r, newParentID)); // Segment goes through current parent
				} else
				{
					segmentsForParents.push_back(Segment(l, r, segmentForChild.newID)); // Segment goes straight from parent of parent to child of parent
				}
				
			}


			/////////////////////////
			// Deal with ancestors //
			/////////////////////////

			if (isAncestor.size() > parentID && isAncestor[parentID])
			{
				if (newParentID == -1) // If the ancestor has not been added to the tree
				{
					/*
						propagate mutations from oldID 'parentID' to newIDs stored in 'segmentsForChildren'
					*/

					for (auto& segment : copy_segmentsForChildren)
					{
						assert(segment.newID < No.size());
						assert(ancestorsGenetics.size() > parentID);
						auto it = std::lower_bound(new_ancestorsID.begin(), new_ancestorsID.end(), segment.newID);
						if (it == new_ancestorsID.end())
						{
							// Not yet registered as ancestor
							NodeGenetics newAncestorGenetics(Ni[parentID].birth);
							newAncestorGenetics.propagateMutationsForSegment(ancestorsGenetics[parentID], segment.left, segment.right);

							new_ancestorsID.insert(it, segment.newID);
							new_ancestorsGenetics.push_back(newAncestorGenetics);
						} else
						{
							auto new_ancestor_index = it - new_ancestorsID.begin();
							assert(new_ancestor_index < new_ancestorsGenetics.size());
							new_ancestorsGenetics[new_ancestor_index].propagateMutationsForSegment(ancestorsGenetics[parentID], segment.left, segment.right);
						}
					}
				}
			}
			segmentsForChildren.resize(0); segmentsForChildren.shrink_to_fit();

			// Reset the list of segments to give to parents
			simplify_mergeSegmentsGoingThroughNode(segmentsForParents);
			listSegments[parentID] = segmentsForParents;
		}
	

		// uppdate childGeneration
		if (edge_index >= 0)
		{
			assert(childGeneration > Ni[Ni_indexForParents].birth);
			childGeneration = Ni[Ni_indexForParents].birth;
		}
	}

	

	/////////////////////////
	// Recreate isAncestor //
	/////////////////////////

	auto max = *std::max_element(new_ancestorsID.begin(), new_ancestorsID.end());
	isAncestor.clear();
	isAncestor.resize(max, false);
	isAncestor.shrink_to_fit();
	for (auto& id : new_ancestorsID)
	{
		isAncestor[id] = true;
	}


	/////////////////////////
	// Swap ancestors info //
	/////////////////////////

	new_ancestorsGenetics.swap(ancestorsGenetics);
	new_ancestorsID.swap(ancestorsID);


	//////////////////////////////////////////////
	// Reverse nodes witout messing up with IDs //
	//////////////////////////////////////////////


	std::cout << "I think it actually messes up the IDs\n";
	No.reverseData();
}


void drawSegments(std::vector<Segment>& segments)
{
	for (auto it = segments.cbegin(); it < segments.cend(); ++it)
	{
		auto l = it->left;
		auto r = it->right;
		auto ID = it->newID;

		assert(l < r);

		std::cout << ID << ":";

		size_t i = 0;
		for ( ; i < l; ++i)
			std::cout << " ";
		for (; i < r; ++i)
			std::cout << "-";
		std::cout << "\n";
	}
}


void T4TreeRec::print()
{
	nodes.print();
	edges.print();
}































/*
	std::vector<std::vector<size_t>> childrenPerHaplos; // childrenPerHaplos[parent_haploID][child_index] = child_haploID



	////////////////////////////////////////
	// Going up the tree to list children //
	////////////////////////////////////////
	{

		// Initialize queue
		std::queue<int> Q;
		size_t lastGeneration = E.back().generation;
		for (int haploID = E.size() - 1 ; haploID >= 0 && E[haploID].generation == lastGeneration ; --haploID)
		{
			Q.push(haploID);
		}

		// Resize childrenPerHaplos
		childrenPerHaplos.resize(E.size() - Q.size());

		// look up the tree
		while (Q.size())
		{
			auto haplo_ID = Q.pop();

			auto parent1 = E[haploID].parents.first;
			if (parent1 != -1)
			{
				childrenPerHaplos[parent1].push_back(haploID);
				Q.push(parent1);
			}
			if (E[haploID].RecPoints.size())
			{
				auto parent2 = E[haploID].parents.second;
				if (parent2 != -1)
				{
					childrenPerHaplos[parent2].push_back(haploID);
					Q.push(parent2);
				}
			}
		}
	}
	

	//////////////////////////////////////////////
	// Going down the tree to build new entries //
	//////////////////////////////////////////////	
	{
		class QElement
		{
		public:
			int parent;
			std::vector<int> children;
			//QElement(int p, std::vector<int>& c):parent(p),children(c){}
		};



		// Initialize queue
		std::queue<QElement> Q;
		for (int haploID = 0 ; haploID < E.size() && E[haploID].parents.first == -1 ; ++haploID)
		{
			Q.push({haploID, childrenPerHaplos[haplo_ID]});
		}

		// Look down the tree and create new entries
		while (Q.size())
		{
			auto Qelement = Q.pop();
			assert(Qelement.children.size());
			
			if (Qelement.children.size() == 1)
			{
				Q.push({
					Qelement.parent,
					childrenPerHaplos[Qelement.children[0]] // grand children
				});
			} else
			{
				for (auto& child_ID : Qelement.children)
				{
					E[child_ID]
					...
				}
			}
		}
	}

	














	std::vector<T4TreeRecEntry> E;
	E.swap(entries);


	std::vector<std::vector<size_t>> childrenPerHaplos(E.size()); // I actually don't need to initialize the last generation. I could use SSP for that. childrenPerHaplos[parent_haploID][child_index] = child_haploID
	std::vector<size_t> whoIsParent(100);

	size_t lastGeneration = E.back().generation;
	for (size_t child_haploID = E.size() - 1 ; child_haploID >= 0 && E[child_haploID].generation == lastGeneration ; --child_haploID)
	{
		auto& firstParentID = E[child_haploID].parents.first;
		childrenPerHaplos[firstParentID].push_back(child_haploID);
		whoIsParent.push_back(firstParentID);
		if (E[haploID].RecPoints.size())
		{
			auto& secondParentID = E[child_haploID].parents.second;
			childrenPerHaplos[secondParentID].push_back(child_haploID);
			whoIsParent.push_back(secondParentID);
		}
	}





	while (whoIsParent.size())
	{
		std::vector<size_t> whoIsParent_b;
		whoIsParent_b.reserve(whoIsParent.size());
		for (auto& parent : whoIsParent)
		{
			if (childrenPerHaplos[parent].size() == 1)
			{
				auto& child_haploID = childrenPerHaplos[parent][0];
				auto& firstParentID = E[child_haploID].parents.first;
				whoIsParent_b.push_back(firstParentID);
				childrenPerHaplos[firstParentID].push_back(child_haploID);
				if (E[haploID].RecPoints.size())
				{
					auto& secondParentID = E[child_haploID].parents.first;
					whoIsParent_b.push_back(secondParentID);
					childrenPerHaplos[secondParentID].push_back(child_haploID);
				}
				childrenPerHaplos[parent].resize(0);
			} else
			{
				assert(childrenPerHaplos[parent].size() > 1);

				for (auto& child_haploID : childrenPerHaplos[parent])
				{

				}

				// Assume overlap rather than giguring it out
				
				// Are they actually all independent as they inheritted different bits of DNA?
				bool isOverlap = false;
				for (size_t child_index = 1 ; child_index < childrenPerHaplos[parent].size() ; ++child_index)
				{
					auto& child_haploID = childrenPerHaplos[parent][child_index];
					if (E[child_haploID].RecPoints.size() == 0)
					{
						isOverlap = true;
						goto OverlapFiguredOut;
					}
					bool whichIsParent_child = E[child_haploID].parents.first == parent ? false : true;
					for (size_t prev_child_index = 0 ; prev_child_index < child_index ; ++prev_child_index)
					{
						auto& prev_child_haploID = childrenPerHaplos[parent][prev_child_index];
						bool whichIsParent_prev_child = E[prev_child_haploID].parents.first == parent ? false : true;
						assert(E[prev_child_haploID].parents[whichIsParent_prev_child] == E[child_haploID].parents[whichIsParent_child]);
						if ()

					}
				}
				OverlapFiguredOut:
				
			}
		}
		whoIsParent.swap(whoIsParent_b);
	}










	// Create stack with the last generation
	std::queue<size_t> Q;
	auto lastGeneration = E.back().generation;
	for (size_t ID = E.size() - 1 ; E[ID].generation == lastGeneration ; --ID)
	{
		Q.push(ID);
	}
	

	while (Q.size())
	{
		auto haploID = Q.top();
		auto& entry = E[haploID];
		auto nbChildren = entry.children.size();

		if (nbChildren == 1)
		{
			auto& child = entry.children.back();
			if (child.RecPoints.size() == 0 )
			{
				// everything comes from one first grandparent
				T4Entry(
					child.generation,
					child.RecPoints,
					std::pair<
						entry.parents.first,
						0 // Would love to put NAN but size_t does not have a NAN value
					>,
					child.children
				);
			}
			(void) Q.pop();
		} else
		{

		}
	}*/


