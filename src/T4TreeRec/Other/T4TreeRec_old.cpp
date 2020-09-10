


void T4TreeRec::initialize(Pop& pop)
{
	/////////////////
	// Reserve RAM //
	/////////////////

	{
		// How much RAM to reserve?
		unsigned long long maxEverNbNodesProduced = 0;
		assert(GP->__GenerationChange.size() == SSP->__patchCapacity.size());
		for (uint32_t generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; ++generation_index )
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
		nodes.data.reserve(nbNodesToReserveFor);
		edges.data.reserve(nbNodesToReserveFor * (1 + SSP->TotalRecombinationRate));
	}
	

/*	/////////////////////////////////
	// Add Ancestor and set pop ID //
	/////////////////////////////////

	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber; ++patch_index )
	{
		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index )
		{
			for (uint32_t haplo_index = 0 ; haplo_index < 2; ++haplo_index )
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

void T4TreeRec::setPopToUniqueID(Pop& pop, uint32_t T4ID)
{
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber; ++patch_index )
	{
		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index )	
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
	std::vector<uint32_t> nbNodesPerGeneration(GP->CurrentGeneration+1, 0);
	for (uint32_t u = 0 ; u < nodes.size() ; ++u )
	{
		assert(nodes[u].birth >= 0 && nodes[u].birth < nbNodesPerGeneration.size());
		++nbNodesPerGeneration[nodes[u].birth];
	}

	uint32_t sum = 0;
	for (uint32_t gen = 0 ; gen < nbNodesPerGeneration.size() ; ++gen)
	{
		std::cout << "Nb nodes in generation: " <<gen << " is " << nbNodesPerGeneration[gen] << "\n";
		sum += nbNodesPerGeneration[gen];
	}
	assert(sum == nodes.size());


	// What fraction of haplotypes have parented at each generation
	std::unordered_map<uint32_t, std::vector<uint32_t>> hash;
	for (uint32_t edge_index = 0 ; edge_index < edges.size() ; ++edge_index)
	{
		hash[edges[edge_index].parent].push_back(edge_index);
	}

	std::vector<uint32_t> nbHaplosThatReproducedAtEachGeneration(GP->CurrentGeneration+1, 0);
	for (uint32_t u = 0 ; u < nodes.size() ; ++u )
	{
		if (hash[u].size())
		{
			++nbHaplosThatReproducedAtEachGeneration[nodes[u].birth];
		}
	}

	assert(nbHaplosThatReproducedAtEachGeneration.size() == nbNodesPerGeneration.size());
	for (uint32_t gen = 0 ; gen < nbNodesPerGeneration.size() ; ++gen)
	{
		assert(nbHaplosThatReproducedAtEachGeneration.size() > gen);
		assert(nbNodesPerGeneration.size() > gen);
		//double frac = (double)nbHaplosThatReproducedAtEachGeneration[gen] / (double)nbNodesPerGeneration[gen];

		std::cout << "Fractions of haplodes that reproduced at generation: " <<gen << " is " << (double)nbHaplosThatReproducedAtEachGeneration[gen] / (double)nbNodesPerGeneration[gen] << "\n";
	}

	

	// Pop IDs
	std::vector<uint32_t> IDtable;
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
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

	/*for (uint32_t ID = 0 ; ID < IDtable.size() ; ++ID)
	{
		std::cout << "popID: " << ID << ": " << IDtable[ID] << "\n";
	}	*/

}


/*class mapT4TreeStructure
{
	std::vector<uint32_t> birthTimeBoundaries;

	void compute(NodeTable& nodes)
	{
		birthTimeBoundaries.resize(0);

		for (uint32_t nodeID = 0 ; nodeID < nodes.size() ; ++nodeID)
		{

		}
	}

	uint32_t beginGeneration()
	{

	}

	uint32_t endGeneration()
	{

	}


};
*/


std::vector<std::vector<std::vector<uint32_t>>> T4TreeRec::placeMutations(Pop& pop, bool isNeedSimplify, bool gatherAndOutputData) // Pop is used to access haplotypes T4ID to do a good matching and to reset them
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

		std::vector<std::vector<std::vector<uint32_t>>> ret(GP->PatchNumber);
		for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
		{
			auto& patch = pop.getPatch(patch_index);
			ret[patch_index].reserve(2*SSP->patchSize[patch_index]);

			assert(patch.getpatchCapacity() >= SSP->patchSize[patch_index]);
			for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
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



	assert(SSP->T4_MutationRate.size() == 1 || SSP->T4_MutationRate.size() == SSP->Gmap.T4_nbLoci);

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

	std::unordered_map<uint32_t, std::vector<uint32_t>> hash;
	for (uint32_t edge_index = 0 ; edge_index < edges.size() ; ++edge_index)
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
	for (uint32_t ancestor_index = 0 ; ancestor_index < ancestorsID.size() ; ++ancestor_index)
	{
		assert(treeGenetics.size() > ancestorsID[ancestor_index]);
		treeGenetics[ancestorsID[ancestor_index]] = ancestorsGenetics[ancestor_index];
		//treeGenetics[ancestorsID[ancestor_index]].whatHasBeenSetYet.push({0, (uint32_t) SSP->Gmap.T4_nbLoci,-1});
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

	uint32_t nbIndsBornInCurrentGeneration = 0; // Just a security
	std::vector<uint32_t> oldToNewIDMap(oldNodes.size());
	for (uint32_t u = 0 ; u < oldNodes.size() ; ++u )
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
			uint32_t v;
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


	std::vector<std::vector<std::vector<uint32_t>>> ret;
	if (gatherAndOutputData)
	{
		////////////////////////////////////////////////////////
		//// Gather data in convenient format and reset IDs ////
		////////////////////////////////////////////////////////

		//std::cout << "About to gather data\n";

		//std::cout << "Finish gathering data\n";

		ret.resize(GP->PatchNumber);
		for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
		{
			auto& patch = pop.getPatch(patch_index);
			ret[patch_index].reserve(2*SSP->patchSize[patch_index]);

			for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
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
	}


	return ret;
}


int T4TreeRec::searchForClonesAmongAncestors(const NodeGenetics& input) const
{
	assert(ancestorsGenetics.size() == ancestorsID.size());

	for (uint32_t ancestor_index = 0 ; ancestor_index < ancestorsGenetics.size() ; ++ancestor_index)
	{
		const auto& ancestor = ancestorsGenetics[ancestor_index];
		
		// test if ancestor is clone of input
		if (input.mutations.size() != ancestor.mutations.size())
		{
			continue;
		} else
		{
			bool areClones = true;
			for (uint32_t mut_index = 0 ; mut_index < input.mutations.size() ; ++mut_index )
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
	assert(nextLeft == SSP->Gmap.T4_nbLoci);
}
*/

void NodeGenetics::propagateMutationsForSegment(const NodeGenetics& parent, const uint32_t left, const uint32_t right)
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
			uint32_t prev = std::numeric_limits<uint32_t>::max();
			for (auto& m : mutations)
			{
				if (prev != std::numeric_limits<uint32_t>::max() && prev >= m)
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

uint32_t T4TreeRec::addAncestor(NodeGenetics& newAncestor)
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

uint32_t T4TreeRec::addHaplotype(std::vector<uint32_t>& RP, std::pair<uint32_t, uint32_t> p)
{
	/*
		Haplotypes need to know their nodeID
	*/


	if (p.first == p.second)
	{
		RP.resize(1);
		RP[0] = std::numeric_limits<uint32_t>::max();
	}

	//std::cout << "p.first = " << p.first << "\n";
	//std::cout << "nodes.size() = " << nodes.size() << "\n";
	assert(nodes.size() > p.first);
	assert(RP.back() == std::numeric_limits<uint32_t>::max());


	
	auto offNodeID = nodes.addNode(GP->CurrentGeneration);
	
	if (RP.size() > 1) // This if else is just for performance
	{
		uint32_t l = 0;
		bool isPfirst = true;
		for (uint32_t DNAsegment = 0 ; DNAsegment < RP.size() ; ++DNAsegment)
		{
			uint32_t r = RP[DNAsegment];
			if (r == std::numeric_limits<uint32_t>::max())
			{
				r = SSP->Gmap.T4_nbLoci;
			} else
			{
				r = SSP->Gmap.FromLocusToNextT4Locus(r); 
			}
			edges.addEdge(l, r, isPfirst ? p.first : p.second, offNodeID);
			l = r;
			isPfirst = !isPfirst;
		}
	} else
	{
		edges.addEdge(0, SSP->Gmap.T4_nbLoci, p.first, offNodeID);
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

void T4TreeRec::simplify(Pop& pop)
{
	lastGenerationSimplified = GP->CurrentGeneration;
	auto& Ei = edges;
	auto& Ni = nodes;
	NodeTable No;
	EdgeTable Eo;

	///////////////
	// Invert Ni //
	///////////////

	//std::reverse(Ni.data.begin(), Ni.data.end());

	////////////////////////
	// Prepare isAncestor //
	////////////////////////

	if (isAncestor.size() < Ni.size())
	{
		isAncestor.resize(Ni.size(), false);
	}



	///////////////////////////////////////////
	// Hash table for finding parental edges //
	///////////////////////////////////////////

	//std::cout << "Hash table\n";
	std::unordered_map<uint32_t, std::vector<uint32_t>> hash;
	for (uint32_t edge_index = 0 ; edge_index < Ei.size() ; ++edge_index)
	{
		hash[Ei[edge_index].parent].push_back(edge_index);
	}


	//////////////
	// Do stuff //
	//////////////

	// Initialize queue
	std::priority_queue<Segment, std::vector<Segment>, std::greater<Segment>> Q;

	// Initialize A (A stands for Ancestral) and map of IDs
	std::vector<std::vector<Segment>> A(Ni.size());
	std::vector<uint32_t> oldToNewIDMap(Ni.size(), std::numeric_limits<uint32_t>::max()); // oldToNewIDMap[oldID] = newID

	//nodes.print();

	uint32_t lastGeneration = Ni.back().birth;
	assert(lastGeneration == GP->CurrentGeneration);
	{
		uint32_t nbNodesInCurrentGeneration = 0;
		for (uint32_t u = Ni.size() - 1 ; u < Ni.size() && Ni[u].birth == lastGeneration ; --u)
		{
			auto v = No.addNode(Ni[u]);
			oldToNewIDMap[u] = v;
			//std::cout << "setting oldToNewIDMap["<<Ni.size()-1-u<<"] = " << v << "\n";
			A[u].push_back(Segment(0, SSP->Gmap.T4_nbLoci, v));
			++nbNodesInCurrentGeneration;
		}
		//std::cout << "nbNodesInCurrentGeneration = " << nbNodesInCurrentGeneration << "\n";
		//std::cout << "SSP->TotalpatchSize = " << SSP->TotalpatchSize << "\n";
		assert(nbNodesInCurrentGeneration == 2*SSP->TotalpatchSize);
	}
	

	// Loop through nodes
	assert(Ni.size());
	for (uint32_t u = Ni.size() - 1 ; u < Ni.size() ; --u)
	{
		// Select parental edges
		auto& parentEdgeIndices = hash[u]; // The naming is poor. 'parentEdgeIndices' is actually a vector to the edge indices where u is parent
		int v = -1;

		// Insert edge ancestry intersections
		for (auto& edge_index : parentEdgeIndices)
		{
			assert(Ei.size() > edge_index);
			auto edge = Ei[edge_index];
 
			assert(A.size() > edge.child);
			for (auto& x : A[edge.child]) // x is therefore a Segment
			{
				if (x.right > edge.left && edge.right > x.left) // intersection
				{
					auto maxleft = std::max(x.left, edge.left);
					auto minright = std::min(x.right, edge.right);
					assert(maxleft <= minright);
					Q.push(Segment(maxleft, minright, x.node));
				}
			}
		}

		
		//std::cout << "For u = " << u << ", Q.size() = "<< Q.size() << "\n";
		
		while (Q.size())
		{
			// Find segments with minimum left coordinate
			auto l = Q.top().left;
			uint32_t r = SSP->Gmap.T4_nbLoci;
			std::vector<Segment> X;
			while (Q.size() && Q.top().left == l)
			{
				auto x = Q.top(); Q.pop();
				X.push_back(x);
				r = std::min(r, x.right);
			}

			if (Q.size())
			{
				r = std::min(r, Q.top().left);
			}

			Segment alpha;
			if (!(X.size() > 1 || isAncestor[u])) // Original was "if (X.size() == 1)"
			{
				// Only one segment starting at minimal left
				auto& x = X[0];
				alpha = x;
				if (Q.size() && Q.top().left < x.right)
				{
					alpha = Segment(x.left, Q.top().left, x.node);
					x.left = Q.top().left;
					Q.push(x);
				}
			} else
			{
				// Overlap new output node
				if (v == -1)
				{
					//std::cout << "Ni[u].birth = " << Ni[u].birth << "\n";
					assert(Ni.size() > u);
					v = No.addNode(Ni[u]);
					oldToNewIDMap[u] = v;
					//std::cout << "setting oldToNewIDMap["<<u<<"] = " << v << "\n";
				}
				
				// Record edges
				alpha = Segment(l,r,v);
				for (auto& x : X)
				{
					Eo.addEdge(l, r, v, x.node);
					if (x.right > r)
					{
						x.left = r;
						Q.push(x);
					}
				}
			}
				
			// Left coordinate loop
			A[u].push_back(alpha);
		}
	}

	/*
	std::cout << "Ei.size() = " << Ei.size() << "\n";
	std::cout << "Eo.size() = " << Eo.size() << "\n";
	std::cout << "Ni.size() = " << Ni.size() << "\n";
	std::cout << "No.size() = " << No.size() << "\n";
	*/


	////////////////
	// Sort nodes //
	////////////////

	nodes.clear(); // remember Ni is nodes. So I am refilling nodes here.
	for (int newID = No.size()-1 ; newID >= 0 ; --newID)
	{
		(void) nodes.addNode(No[newID]); // returns newnewID (which is equal to nodes.size() - 1 - newID) but I don't need it here
	}
	assert(No.size() == nodes.size());

	// Do not use Ni as Ni now contains the sorted No and the name Ni would be very confusing.


	///////////////////
	// Compact edges //
	///////////////////

	std::sort(
		Eo.data.begin(),
		Eo.data.end(), 
		std::greater<Edge>()
	);
	
	edges.clear(); // Remember Ei is edges so I am filling up edges here
	uint32_t start = 0;
	for (uint32_t j = 1 ; j < Eo.size(); ++j)
	{
		bool condition =
			Eo[j-1].right != Eo[j].left ||
			Eo[j-1].parent != Eo[j].parent ||
			Eo[j-1].child != Eo[j].child
			;
		if (condition)
		{
			edges.addEdge(
				Eo[start].left,
				Eo[j-1].right,
				nodes.size() - 1 - Eo[j-1].parent,
				nodes.size() - 1 - Eo[j-1].child
			);
			start = j;
		}
	}
	if (Eo.size())
	{
		edges.addEdge(
			Eo[start].left,
			Eo.back().right,
			nodes.size() - 1 - Eo.back().parent,
			nodes.size() - 1 - Eo.back().child
		);
	}

	// Do not use Ei as Ei now contains the sorted and compacted Eo and the name Ei would be very confusing.

	///////////////////////////
	// Redefine ancestors ID //
	///////////////////////////

	assert(ancestorsID.size());
	assert(ancestorsID.size() == ancestorsGenetics.size());

	std::vector<uint32_t> new_ancestorsID;
	std::vector<NodeGenetics> new_ancestorsGenetics;

	uint32_t prevAncestorID = std::numeric_limits<uint32_t>::max() ; // just a security
	for (uint32_t ancestor_index = 0 ; ancestor_index < ancestorsID.size() ; ++ancestor_index)
	{
		//std::cout << "ancestorsID["<<ancestor_index<<"] = " << ancestorsID[ancestor_index] << "\n";
		auto& u = ancestorsID[ancestor_index];
		assert(prevAncestorID != u); // just a security. Must be true because I don't keep clones
		if (oldToNewIDMap[u] != std::numeric_limits<uint32_t>::max())
		{
			assert(oldToNewIDMap[u] < nodes.size());
			auto newnewv = nodes.size() - 1 - oldToNewIDMap[u];
			prevAncestorID = u; // just a security
			new_ancestorsID.push_back(newnewv);
			new_ancestorsGenetics.push_back(ancestorsGenetics[ancestor_index]);
		}
	}
	assert(new_ancestorsID.size());
	assert(new_ancestorsID.size() == new_ancestorsGenetics.size());
	ancestorsID.swap(new_ancestorsID);
	ancestorsGenetics.swap(new_ancestorsGenetics);


	//////////////////////////
	// Reset T4flags in Pop //
	//////////////////////////

	for (uint32_t patch_index = 0 ; patch_index < pop.getNbPatches() ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);
		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			auto& ind = patch.getInd(ind_index);

			auto& haplo0 = ind.getHaplo(0);
			auto& haplo1 = ind.getHaplo(1);

			/*
			std::cout << "haplo0.T4ID = "<< haplo0.T4ID <<"\n";
			std::cout << "haplo1.T4ID = "<< haplo1.T4ID <<"\n";
			std::cout << "oldToNewIDMap.size() = "<< oldToNewIDMap.size() <<"\n";
			*/

			assert(oldToNewIDMap.size() > haplo0.T4ID);
			assert(oldToNewIDMap.size() > haplo1.T4ID);

			auto h0newID = oldToNewIDMap[haplo0.T4ID];
			auto h1newID = oldToNewIDMap[haplo1.T4ID];
			//std::cout << "Withdrawing newID = "<< h0newID << " from oldToNewIDMap with oldID "<< haplo0.T4ID <<"\n";
			//std::cout << "Withdrawing newID = "<< h1newID << " from oldToNewIDMap with oldID "<< haplo1.T4ID <<"\n";

			auto h0newnewID = nodes.size() - h0newID - 1; // That mess is because nodes have been reverse and therefore new is not new anymore!
			auto h1newnewID = nodes.size() - h1newID - 1; // That mess is because nodes have been reverse and therefore new is not new anymore!

			// Make sure it did not overflow below zero (it is an unsigned type)
			assert(h0newnewID < nodes.size());
			assert(h1newnewID < nodes.size());


			haplo0.T4ID = h0newnewID;
			haplo1.T4ID = h1newnewID;
		}	
	}
	

	////////////////////////////////////////////////////
	// Place mutations if too many edges to ancestors //
	////////////////////////////////////////////////////

	//if (nodes.size() > 50 &&  && (nodes.back().birth - nodes[0].birth > 100))
}


void T4TreeRec::print()
{
	nodes.print();
	edges.print();
}

std::vector<uint32_t> T4TreeRec::getMutationsOfID(const uint32_t id) const
{
	auto i = std::find(ancestorsID.begin(), ancestorsID.end(), id) - ancestorsID.begin();
	assert(i >= 0 && i < ancestorsGenetics.size());
	return ancestorsGenetics[i].mutations;
}



























/*
	std::vector<std::vector<uint32_t>> childrenPerHaplos; // childrenPerHaplos[parent_haploID][child_index] = child_haploID



	////////////////////////////////////////
	// Going up the tree to list children //
	////////////////////////////////////////
	{

		// Initialize queue
		std::queue<int> Q;
		uint32_t lastGeneration = E.back().generation;
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


	std::vector<std::vector<uint32_t>> childrenPerHaplos(E.size()); // I actually don't need to initialize the last generation. I could use SSP for that. childrenPerHaplos[parent_haploID][child_index] = child_haploID
	std::vector<uint32_t> whoIsParent(100);

	uint32_t lastGeneration = E.back().generation;
	for (uint32_t child_haploID = E.size() - 1 ; child_haploID >= 0 && E[child_haploID].generation == lastGeneration ; --child_haploID)
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
		std::vector<uint32_t> whoIsParent_b;
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
				for (uint32_t child_index = 1 ; child_index < childrenPerHaplos[parent].size() ; ++child_index)
				{
					auto& child_haploID = childrenPerHaplos[parent][child_index];
					if (E[child_haploID].RecPoints.size() == 0)
					{
						isOverlap = true;
						goto OverlapFiguredOut;
					}
					bool whichIsParent_child = E[child_haploID].parents.first == parent ? false : true;
					for (uint32_t prev_child_index = 0 ; prev_child_index < child_index ; ++prev_child_index)
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
	std::queue<uint32_t> Q;
	auto lastGeneration = E.back().generation;
	for (uint32_t ID = E.size() - 1 ; E[ID].generation == lastGeneration ; --ID)
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
						0 // Would love to put NAN but uint32_t does not have a NAN value
					>,
					child.children
				);
			}
			(void) Q.pop();
		} else
		{

		}
	}*/


