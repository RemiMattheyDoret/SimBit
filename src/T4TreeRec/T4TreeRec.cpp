
void T4TreeRec::shift_generations_after_burn_in()
{
	int maxGen = std::numeric_limits<int>::lowest();
	//std::cout << "maxGen = " << maxGen << "\n";
	for (uint32_t nodeID = 0 ; nodeID < nodes.size() ; ++nodeID)
	{
		//std::cout << "bb " << nodes[nodeID].generation << "\n";
		maxGen = std::max(maxGen, nodes[nodeID].generation);
		//std::cout << "maxGen = " << maxGen << "\n";
	}
	if (maxGen >= 0)
	{
		std::cout << "It is likely an internal error but it looks like the neutral burn in lasted more than (the absolute of) " << std::numeric_limits<int>::lowest() << " generations. maxGen = "<< maxGen <<". The burn in cannot be any longer!\n";
		abort();
	}


	for (uint32_t nodeID = 0 ; nodeID < nodes.size() ; ++nodeID)
	{
		assert(nodes[nodeID].generation < 0);
		//std::cout << nodes[nodeID].generation << " -> " << nodes[nodeID].generation - maxGen + GP->startAtGeneration <<  "\n";
		nodes[nodeID].generation = nodes[nodeID].generation - maxGen + GP->startAtGeneration;		
	}
	assert(nodes.back().generation == GP->startAtGeneration);
	ancestralGeneration = ancestralGeneration - maxGen + GP->startAtGeneration;
	assert(ancestralGeneration < 0);

	for (NodeGenetics& ancNode : ancestorsGenetics)
	{
		ancNode.setGeneration(ancestralGeneration);
	}
}


void T4TreeRec::initialize(Pop& pop)
{
	haveAllLociCoallesced_info = 'u';
	/////////////////
	// Reserve RAM //
	/////////////////

	{
		// How much RAM to reserve?
		

		int totalInitialPatchCapacity = 0;
		for (auto& pc : SSP->__patchCapacity[0])
		{
			totalInitialPatchCapacity += pc;
		}


		double nbNodes = 2.0 * (double) SSP->T4_simplifyEveryNGenerations * (double)totalInitialPatchCapacity;
		nbNodes = nbNodes > 5e8 ? 5e8 : nbNodes;
		double nbEdges = 2.0 * (double) nbNodes * (double)(1.0 + SSP->TotalRecombinationRate);
		nbEdges = nbEdges > 1e9 ? 1e9 : nbEdges;

		{
			auto& x = nbNodes;
			x += 10 * pow(x, 0.5);
			x += 2 * x / SSP->T4_simplifyEveryNGenerations;
		}
		{
			auto& x = nbEdges;
			x += 10 * pow(x, 0.5);
			x += 2 * x / SSP->T4_simplifyEveryNGenerations;
		}


		// Reserve RAM		
		//std::cout << "reserving " << nbEdges << " edges\n";
		edges.reserve(nbEdges);
		
		//std::cout << "reserving " << nbNodes << " nodes\n";
		nodes.reserve(nbNodes ); // Add one generation worth of nodes

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


	///////////////////
	// Add Ancestors //
	///////////////////

	ancestralGeneration = GP->burnInUntilT4Coal_check_every_N_generations == -1 ? GP->startAtGeneration : std::numeric_limits<int>::lowest() ;
	NodeGenetics ancestorGenetics; // initialize at 0 for all loci

	for (uint32_t patch_index = 0 ;  patch_index < GP->PatchNumber; ++patch_index)
		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index)
			for (uint32_t haplo_index = 0 ; haplo_index < 2; ++haplo_index)
				pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).T4ID = addAncestor(ancestorGenetics, patch_index, GP->burnInUntilT4Coal_check_every_N_generations == -1 ? GP->startAtGeneration : std::numeric_limits<int>::lowest());
	
	/*	
		auto T4ID = addAncestor(ancestorGenetics); // Just one node is enough as they are all clones
		setPopToUniqueID(pop, T4ID);
	*/
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


std::vector<std::vector<std::vector<uint32_t>>> T4TreeRec::getMutations(Pop& pop)
{
	assert(nodes.size() == ancestorsGenetics.size());

	assert(isAlreadyAsSimplifiedAsPossible);

	assert(nodes.size() == ancestorsGenetics.size());

	std::vector<std::vector<std::vector<uint32_t>>> ret(GP->PatchNumber);
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);
		ret[patch_index].reserve(2*SSP->patchSize[patch_index]);

		assert(patch.getpatchCapacity() >= SSP->patchSize[patch_index]);
		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
			{
				auto T4ID = patch.getInd(ind_index).getHaplo(haplo_index).T4ID;

				// Build ret (object to return)
				ret[patch_index].push_back(getMutationsOfID(T4ID));
			}
		}	
	}


	return ret;
}


/*
void T4TreeRec::setPartOfNewAncestor(std::vector<NodeGenetics>& allNodes, std::vector<NodeGenetics>& newAncestorsGenetics)
{
	assert(newAncestorsGenetics.size() == SSP->TotalpatchSize * 2);
	for (int newID = 0 ; newID < SSP->TotalpatchSize * 2 ; ++newID)
	{
		auto oldID = allNodes.size() - newID - 1;
		newAncestorsGenetics[newID].push_backMutations(allNodes[oldID].moveMutations());
	}
}
*/

/*
void T4TreeRec::resetAncestorsEdgesAndNodesAfterPlacingMutations(std::vector<NodeGenetics>& allNodes, bool shouldCompleteDeleteAllTree)
{
	NodeTable oldNodes;
	EdgeTable oldEdges;
	oldNodes.swap(nodes);
	oldEdges.swap(edges);
	nodes.reserve(SSP->TotalpatchSize * 2);

	ancestorsGenetics.clear();
	ancestorsGenetics.reserve(SSP->TotalpatchSize * 2);

	assert(nodes.size() == 0);
	assert(edges.size() == 0);
	for (int newID = 0 ; newID < SSP->TotalpatchSize * 2 ; ++newID)
	{
		auto oldID = allNodes.size() - newID - 1;
		assert(newID == nodes.addNode(oldNodes[oldID]));
		nodes.back().generation = oldNodes[oldID].generation;
		
		//std::cout << "oldID = " << oldID << "\n";
		//std::cout << "newID = " << newID << "\n";
		//std::cout << "allNodes.size() = " << allNodes.size()  << "\n";
		
		assert(allNodes.size() > oldID);
		allNodes[oldID].resetID(newID);
		if (!shouldCompleteDeleteAllTree)
		{
			ancestorsGenetics.push_back(allNodes[oldID]);
			ancestorsGenetics.back().setGeneration(GP->CurrentGeneration);
		}
	}
	assert(nodes.size() == SSP->TotalpatchSize * 2);
	assert(edges.size() == 0);
	ancestralGeneration = GP->CurrentGeneration;
}
*/


void T4TreeRec::assertIsFullySimplified() const
{
	assert(isAlreadyAsSimplifiedAsPossible);
}


std::pair< std::vector<std::vector<std::vector<uint32_t>>>, std::vector<T4TreeRec::PatchCountData> > T4TreeRec::placeMutations(Pop& pop, bool shouldDeleteTree, bool shouldCompleteDeleteAllTree, unsigned char outputType) // Pop is used to access haplotypes T4ID to do a good matching and to reset them
{
/*
	It computes the current states and set all current individuals as ancestors (unless if there are clones, then they share the same T4ID)
*/
	/*
	std::cout << "Entering place mutations " << nodes.size() << " nodes and " << edges.size() << " edges\nEntering tree:\n";
	this->print();
	std::cout << "###############################\n\n\n";
	*/
	if (shouldDeleteTree || shouldCompleteDeleteAllTree) assert(shouldDeleteTree != shouldCompleteDeleteAllTree);
	assert(outputType == 'r' || outputType == 'f' || outputType == 'b'); // Raw data, Frequency data, or Both


	std::vector<T4TreeRec::PatchCountData> ret_freqs;
	if (outputType != 'r')
	{
		ret_freqs.resize(GP->PatchNumber, {SSP->T4SNPfreq_assumeMuchDiversity, SSP->Gmap.T4_nbLoci});
	}

	//////////////////////////////////////////////////////////////////////////////
	//// If I don't need to place mutations but just return the current state ////
	//////////////////////////////////////////////////////////////////////////////

	/*
	isAlreadyAsSimplifiedAsPossible = true;
	if (nodes.size() == ancestorsGenetics.size())
	{
		return getMutations(pop);
	}
	*/
	assert(nodes.size());
	assert(SSP->T4_MutationRate.size() == 1 || SSP->T4_MutationRate.size() == SSP->Gmap.T4_nbLoci);


	///////////////////////////////
	//// New ancestor genetics ////
	///////////////////////////////

	std::vector<NodeGenetics> newAncestorsGenetics;
	newAncestorsGenetics.resize( SSP->TotalpatchSize * 2 );


	/////////////////////////////////////////////////////////////////////////
	//// Loop over fractions of genome to consider to minimize RAM usage ////
	/////////////////////////////////////////////////////////////////////////

	size_t fractionGenomeSize = std::ceil((double)SSP->Gmap.T4_nbLoci / (double)SSP->T4_nbRunsToPlaceMutations);
	assert(fractionGenomeSize);
	size_t fractionGenomeFrom = 0;
	size_t fractionGenomeTo = 0;
	for (size_t fractionGenomeIndex = 0 ; fractionGenomeIndex < SSP->T4_nbRunsToPlaceMutations ; ++fractionGenomeIndex)
	{

		fractionGenomeFrom = fractionGenomeTo;
		fractionGenomeTo   = fractionGenomeTo + fractionGenomeSize;
		if (fractionGenomeTo > SSP->Gmap.T4_nbLoci)  fractionGenomeTo = SSP->Gmap.T4_nbLoci;
		//std::cout << "fractionGenomeFrom = " << fractionGenomeFrom << "\n";
		//std::cout << "fractionGenomeTo = " << fractionGenomeTo << "\n";
		assert(fractionGenomeFrom < fractionGenomeTo);
		if (fractionGenomeIndex == SSP->T4_nbRunsToPlaceMutations-1) assert(fractionGenomeTo == SSP->Gmap.T4_nbLoci);




		//////////////////////////////////////////////////////////////////////////////////
		//// Figure at which edge index an individual is producing its last offspring ////
		//////////////////////////////////////////////////////////////////////////////////

		std::vector<bool> isLastReproduction;
		isLastReproduction.resize(edges.size());
		{
			std::vector<bool> hasBeenFoundYet(nodes.size(), false);
			assert(!hasBeenFoundYet[0]);
			assert(hasBeenFoundYet.size() == nodes.size());
			//std::cout << "edges.size() = " << edges.size() << "\n";
			for (size_t edge_index = edges.size()-1 ; edge_index <= edges.size()-1 ; --edge_index)
			{
				if ( edges[edge_index].left < fractionGenomeTo && edges[edge_index].right > fractionGenomeFrom) 
				{
					auto p = edges[edge_index].parent;
					assert(p < nodes.size());

					if (hasBeenFoundYet[p])
					{
						isLastReproduction[edge_index] = false;
					} else
					{
						//std::cout << edge_index << " is last reproduction of "<<p<<" \n";
						hasBeenFoundYet[p] = true;
						isLastReproduction[edge_index] = true;
					}
				}					
			}
		}
		assert(isLastReproduction.size() == edges.size());


		//////////////////////////////////
		//// initialize Tree genetics ////
		//////////////////////////////////

		std::vector<NodeGenetics> allNodes; // This might be improved
		allNodes.resize(nodes.size()); 




		//////////////////////////////////////
		//// Take the ancestral variation ////
		//////////////////////////////////////

		//std::cout << "About to take ancestral variation\n";
		for (uint32_t ancestor_index = 0 ; ancestor_index < ancestorsGenetics.size() ; ++ancestor_index)
		{
			assert(ancestorsGenetics[ancestor_index].getID() < allNodes.size() && ancestorsGenetics[ancestor_index].getID() >= 0);
			allNodes[ancestorsGenetics[ancestor_index].getID()] = ancestorsGenetics[ancestor_index].getSubsetMutations(fractionGenomeFrom, fractionGenomeTo);
			//std::cout << "ancestor: allNodes[" << ancestorsGenetics[ancestor_index].getID()<<"].getGeneration() = " << allNodes[ancestorsGenetics[ancestor_index].getID()].getGeneration() << "\n";
		}
		



		///////////////////////////////////////////
		//// Propagate mutations down the tree ////
		///////////////////////////////////////////

		for (uint32_t edge_index = 0 ; edge_index < edges.size() ; ++edge_index)
		{
			auto& edge = edges[edge_index];
			if ( edges[edge_index].left < fractionGenomeTo && edges[edge_index].right > fractionGenomeFrom) 
			{
				auto& parent = allNodes[edge.parent];

				if (allNodes[edge.child].getID() == -1)
				{
					allNodes[edge.child] = NodeGenetics({}, nodes[edge.child].generation, edge.child);
				}

				
				auto left  = edge.left > fractionGenomeFrom ? edge.left : fractionGenomeFrom;
				auto right = edge.right < fractionGenomeTo ? edge.right : fractionGenomeTo;
				assert (right > left);
					

				propagateMutationsForSegment(allNodes[edge.child], parent, left, right);

				if (isLastReproduction[edge_index])
				{
					parent.clear();
					parent.shrink_to_fit();
				}
			}
		}

		/////////////////////////////////////////////////////////
		//// Compute frequencies if needed and set ancestors ////
		/////////////////////////////////////////////////////////

		if (outputType != 'r') // if not raw output only
		{
			assert(ret_freqs.size() == GP->PatchNumber);
			for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
			{
				auto& patch = pop.getPatch(patch_index);
				for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
				{
					auto& ind = patch.getInd(ind_index);
					for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
					{
						auto oldID = ind.getHaplo(haplo_index).T4ID;
						for (auto& mutatedLocus : allNodes[oldID].getMutations())
						{
							++ret_freqs[patch_index][mutatedLocus];
						}						

						if (outputType != 'f' || shouldDeleteTree)
						{
							auto newID = allNodes.size() - oldID - 1;
							newAncestorsGenetics[newID].push_backMutations(allNodes[oldID].moveMutations());
						}						
					}
				}
			}
		}
	}



		
	////////////////////////////////////////////////////////////
	//// Gather data in convenient format and reset pop IDs ////
	////////////////////////////////////////////////////////////
	
	assert(newAncestorsGenetics.size() == 2 * SSP->TotalpatchSize);
	std::vector<std::vector<std::vector<uint32_t>>> ret;
	if (outputType != 'f') ret.resize(GP->PatchNumber);
		
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);
		if (outputType != 'f') ret[patch_index].reserve(2*SSP->patchSize[patch_index]);

		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			auto& ind = patch.getInd(ind_index);
			for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
			{
				// Get haplotype
				auto& haplo = ind.getHaplo(haplo_index);

				// IDs
				auto oldID = haplo.T4ID;
				assert(oldID < nodes.size());
				auto newID = nodes.size() - oldID - 1;
				assert(newID >= 0 && newID < SSP->TotalpatchSize * 2);
				if (shouldDeleteTree || shouldCompleteDeleteAllTree)
					haplo.T4ID = newID;
				//std::cout << "Set node ID " << newID << ". oldID = " << oldID << "\n";
				newAncestorsGenetics[newID].resetID(newID);
				newAncestorsGenetics[newID].setGeneration(GP->CurrentGeneration);
				

				// return object
				if (outputType != 'f') 
				{
					if (shouldCompleteDeleteAllTree)
					{
						ret[patch_index].push_back(newAncestorsGenetics[newID].moveMutations());
					} else
					{
						ret[patch_index].push_back(newAncestorsGenetics[newID].getMutations());
					}
				}
			}
		}	
	}
	if (outputType != 'f') assert(ret.size() == GP->PatchNumber);

	/////////////////////////
	//// reset ancestors ////
	/////////////////////////
	if (shouldDeleteTree)
	{
		ancestorsGenetics.swap(newAncestorsGenetics);
		ancestorsGenetics.shrink_to_fit();
		std::vector<NodeGenetics>().swap(newAncestorsGenetics);
	} else if (shouldCompleteDeleteAllTree)
	{
		std::vector<NodeGenetics>().swap(ancestorsGenetics);
		std::vector<NodeGenetics>().swap(newAncestorsGenetics);
		this->deleteAll();
	} else
	{
		std::vector<NodeGenetics>().swap(newAncestorsGenetics);
	}

	//std::cout << "Finish gathering data\n";

	//std::cout << "Tree deleted (if asked for)\n";

	/*
	std::cout << "Mutations placed at generation " << GP->CurrentGeneration << ". Ancestral generation is " << ancestralGeneration << " and the tree is\n";
	this->print();
	*/
	//std::cout << "Exit T4TreeRec::placeMutations\n";

	/*
	std::cout << "Exiting place mutations " << nodes.size() << " nodes and " << edges.size() << " edges\nExiting tree:\n";
	this->print();
	std::cout << "###############################\n\n\n";
	*/
	return {ret, ret_freqs};
}


uint32_t T4TreeRec::addAncestor(NodeGenetics& newAncestor, int patch_index, int generation)
{
	auto offNodeID = nodes.addNode({generation, patch_index});
	newAncestor.resetID(offNodeID);
	ancestorsGenetics.push_back(newAncestor);
	return offNodeID;
}


uint32_t T4TreeRec::addHaplotype(std::vector<uint32_t>& RP, std::pair<uint32_t, uint32_t> p, int patch_index)
{
	/*
		Haplotypes need to know their nodeID
	*/

	//if (GP->CurrentGeneration >= 0) std::cout << "Enters T4TreeRec::addHaplotype\n";

	isAlreadyAsSimplifiedAsPossible = false;

	if (p.first == p.second)
	{
		RP.resize(1);
		RP[0] = std::numeric_limits<uint32_t>::max();
	}

	/*
	std::cout << "p.first = " << p.first << "\n";
	std::cout << "p.second = " << p.second << "\n";
	std::cout << "nodes.size() = " << nodes.size() << "\n";
	*/
	if (p.first >= nodes.size() )
	{
		std::cout << "The T4 coalescent tree received a T4ID for a parental haplotype that does not exist in the tree. That might be caused by en internal error or because you set a T4ID for an defining T4IDs yourself such as you can do with --indTypes. T4ID received is " << p.first << " while the current maximal T4ID is " << nodes.size()-1 << "\n";
		abort();
	}
	assert(nodes.size() > p.first);
	assert(RP.back() == std::numeric_limits<uint32_t>::max());


	
	auto offNodeID = nodes.addNode({GP->CurrentGeneration, patch_index});
	assert(p.first < offNodeID);

	if (RP.size() > 1) // This if else is just for performance
	{
		uint32_t l = 0;
		bool isPfirst = true;
		for (uint32_t DNAsegment = 0 ; DNAsegment < RP.size() ; ++DNAsegment)
		{
			uint32_t r;
			if (RP[DNAsegment] == std::numeric_limits<uint32_t>::max())
			{
				r = SSP->Gmap.T4_nbLoci;
			} else
			{
				r = SSP->Gmap.FromLocusToNextT4Locus(RP[DNAsegment]); 
			}

			if (l < r)
			{
				edges.addEdge(l, r, isPfirst ? p.first : p.second, offNodeID);
				l = r;
			}

			if (r == SSP->Gmap.T4_nbLoci)
				break;

			isPfirst = !isPfirst;
		}
	} else
	{
		edges.addEdge(0, SSP->Gmap.T4_nbLoci, p.first, offNodeID);
	}

	//std::cout << "offNodeID = " << offNodeID << "\n";
	//if (GP->CurrentGeneration >= 0) std::cout << "Exit T4TreeRec::addHaplotype\n";

	return offNodeID;
}

void T4TreeRec::simplify_ifNeeded(Pop& pop)
{
	if (
		GP->CurrentGeneration != GP->startAtGeneration &&
		GP->CurrentGeneration % SSP->T4_simplifyEveryNGenerations == 0 &&
		(GP->nbGenerations - GP->CurrentGeneration) > 10
	)
	{
		//std::cout << "Before simplification: edges.size() = "<<edges.size()<<"\n";
		//this->print();
		//auto nbEdgesBefore = edges.size();
		simplify(pop);
		//std::cout << "Simplification went from "<<nbEdgesBefore<<" edges to "<< edges.size()<<" edges.\n";
	}
}

std::vector<T4TreeRec::PaintedSegmentDiversity> T4TreeRec::computeSmallSegmentsDiversityFromPaintedHaplotypes(const std::vector<std::vector<T4TreeRec::PaintedSegment>>& allGatheredSegments) const
{
	/*
	[0 10 0][10 20 1]
	[0 5  0][5  12 2]
	[0 20 0]
	[0 5  0][5  10 0][10 20 1]
	[0 20 3]

	[0  5 F C][5  10 F C][10 15 F C][15 20 F C]

	*/

	/*
	if (GP->CurrentGeneration == 8)
	{
		std::cout << "At generation 8 got the following sequences\n";
		for (const auto& segments : allGatheredSegments)
		{
			for (const auto& segment : segments)
			{
				std::cout << "[" << segment.segment.left << "-" << segment.segment.right << "]\n";
			}
		}
	}
	*/



	// double for faster calculation below but it should be uint32_t
	double nbHaplotypesSampled = allGatheredSegments.size();
	assert(nbHaplotypesSampled > 0.0);

	// Object to return
	std::vector<T4TreeRec::PaintedSegmentDiversity> r;


	// generate all transition points
    std::vector<T4TreeRec::PointForFindingOverlap> points;
    for (const auto& segments : allGatheredSegments)
    {
    	for (const auto& segment : segments)
    	{
    		// Note ind_ID is ancestor_ID
			points.push_back({segment.segment.left, true, segment.ind_ID});
			points.push_back({segment.segment.right, false, segment.ind_ID});
    	}
    }

    // sort transition points
    std::sort(points.begin(), points.end(), 
      [](const T4TreeRec::PointForFindingOverlap& a, const T4TreeRec::PointForFindingOverlap& b) { return a.location < b.location; });


    // initialize overlaps
    std::multiset<int> overs{points[0].ID};


    // for every adjacent transition point
    for (auto i = 1u; i < points.size(); ++i) 
    {
        auto &a = points[i - 1];
        auto &b = points[i];

        // if there is a jump in between transition points
		if (a.location < b.location)
		{
			// Assertions
			assert(overs.size() == nbHaplotypesSampled);
			
			if (r.size())
			{
				assert(r.back().right == a.location);
			} else
			{
				assert(0 == a.location);
			}


			// Add segment ends and diversity (heterozygosity of segments) in returned object
			assert(overs.size());

			/*
			if (GP->CurrentGeneration == 8)
			{
				std::cout << "Between points " << a.location << "-" << b.location << ", overs = [ ";
				for (auto& elem : overs) {std::cout << elem << " ";}
				std::cout << "]\n";
			}
			*/

			if (overs.size() == 1)
			{
				//r.push_back({a.location, b.location, 0.0});
				//r.push_back({a.location, b.location, 1.0, *overs.begin()});
				if (!SSP->SegDiversityFile_includeMainColor && r.size() && r.back().heterozygosity == 0.0)
				{
					assert(r.back().right < b.location);
					r.back().right = b.location;
				} else
				{
					r.push_back({a.location, b.location, 1, *(overs.begin()),  0.0});
				}
			} else
			{
				// Compute counts
				std::vector<std::pair<size_t, int>> counts; // first is count and second is ID
				auto prev = *(overs.begin());
				counts.push_back({0, prev});
				for (auto it = overs.begin() ; it != overs.end() ; ++it)
				{
					if (*it != prev)
					{
						prev = *it;
						counts.push_back({0, prev});
					}
					++(counts.back().first);
				}

				/*
				if (GP->CurrentGeneration == 8)
				{
					for (auto& count : counts)
					{
						std::cout << count.first << ": " << count.second << "\n";
					}
				}
				*/
					

				// Make entries
				/*
				size_t allCounts = 0; // for assertions
				for (auto& count : counts )
				{
					double freq = (double)count.first / (double)nbHaplotypesSampled;
					assert(freq > 0.0 && freq <= 1.0);
					r.push_back({a.location, b.location, freq, count.second});
					allCounts += count.first;
					//std::cout << "count.first = " << count.first << "\n";
					assert(count.first > 0);
				}
				//std::cout << "nbHaplotypesSampled = " << nbHaplotypesSampled << "\n";
				//std::cout << "allCounts = " << allCounts << "\n";
				assert(allCounts == nbHaplotypesSampled);
				*/


				
				// Compute diversity
				size_t sumOfCounts = 0; // For assertion
				double J = 0.0; // J as in Nei's paper (H = 1 - j)
				size_t maxCount = 0; // used only if SSP->SegDiversityFile_includeMainColor
				int mainColor = std::numeric_limits<int>::min();  // used only if SSP->SegDiversityFile_includeMainColor
				for (auto& count : counts)
				{
					sumOfCounts += count.first;
					J += pow(count.first / nbHaplotypesSampled, 2);
					if (maxCount < count.first)
					{
						maxCount = count.first;
						mainColor = count.second;
					}
				}
				assert(maxCount);
				assert(sumOfCounts == nbHaplotypesSampled);
				assert(J >= 0.0 && J <= 1.0);
				auto H = 1-J;
				
				/*
				if (GP->CurrentGeneration == 8)
				{
					std::cout << "H = " << H << "\n";
				}
				*/

				// Add entry in the returned object
				if (!SSP->SegDiversityFile_includeMainColor && r.size() && r.back().heterozygosity == H && r.back().nbColors == counts.size())
				{
					assert(r.back().right < b.location);
					r.back().right = b.location;
				} else
				{

					r.push_back({a.location, b.location, counts.size(), mainColor, H});
				}
			}
		}

        // update overlaps
        if (b.overlap)
           overs.insert(b.ID);
        else
           overs.erase(overs.find(b.ID));  
    }

    assert(r.back().right == SSP->Gmap.T4_nbLoci);
    return r;
}

std::vector<T4TreeRec::PaintedSegmentFrequency> T4TreeRec::computeLargeSegmentsFrequenciesFromPaintedHaplotypes(const std::vector<std::vector<T4TreeRec::PaintedSegment>>& allGatheredSegments) const
{
	// That's a little slow but easy.

	// First linearize the segments
	std::vector<T4TreeRec::PaintedSegment> linearizedPaintedSegments;

	for (uint32_t i = 0 ; i < allGatheredSegments.size() ; ++i)
	{
		linearizedPaintedSegments.reserve(linearizedPaintedSegments.capacity() + allGatheredSegments[i].size());
		for (const auto& elem : allGatheredSegments[i])
			linearizedPaintedSegments.push_back(elem);
	}

	// Then sort them (I could have merge sort them to be faster)
	std::sort(linearizedPaintedSegments.begin(), linearizedPaintedSegments.end());
	if (linearizedPaintedSegments.size() >= 2) assert(linearizedPaintedSegments[0] < linearizedPaintedSegments[1]);
	assert(linearizedPaintedSegments.size());

	// Compute count
	std::vector<T4TreeRec::PaintedSegmentFrequency> data;
	
	auto& prev = linearizedPaintedSegments[0];
	data.push_back({prev.segment.left, prev.segment.right, 0.0});
	for (const auto& elem : linearizedPaintedSegments)
	{
		if (
			elem.segment.left != prev.segment.left ||
			elem.segment.right != prev.segment.right ||
			elem.ind_ID != prev.ind_ID
			)
		{
			data.push_back({elem.segment.left, elem.segment.right, 0.0});
			prev = elem;
		} else assert(elem.patch_index == prev.patch_index);
		++(data.back().freq);
	}

	// Compute relative frequencies
	for (auto& elem : data)
	{
		elem.freq /= allGatheredSegments.size();
		assert(elem.freq >= 0.0 && elem.freq <= 1.0);
	}

	return data;
}



void T4TreeRec::writePaintedHaplotypesDiversity(const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches, const int paintedGeneration, const int observedGeneration, OutputFile& file) const
{
	//std::cout << "enters in T4TreeRec::writePaintedHaplotypesDiversity\n";
	assert(paintedGeneration < observedGeneration);
	assert(focalT4IDs_patches.size() == focalT4IDs.size());

	if (focalT4IDs.size() == 0) return;



	auto allGatheredSegments = computePaintedHaplotypes(paintedGeneration, focalT4IDs, focalT4IDs_patches);
	
	//for (const auto& focalT4ID : focalT4IDs)
	//	allGatheredSegments.push_back(computePaintedHaplotypes(paintedGeneration, focalT4ID));

	///////////////////////////////////////////////////////////
	// Compute frequencies and format into string for return //
	///////////////////////////////////////////////////////////
	/*
		Format: [left right freq][left right freq][left right freq]
		ID is T4ID -> Watch out, one can compare T4IDs of different observedGeneration!
	*/

	assert(focalT4IDs_patches.size());
	int focalT4IDs_patches_index = -1;

	file.open();
	
	//std::cout << "allGatheredSegments.size() = " << allGatheredSegments.size() << "\n";

	for (size_t index_for_patch = 0 ; index_for_patch < allGatheredSegments.size() ; ++index_for_patch)
	{
		// Get patch 
		if (focalT4IDs_patches_index == -1)
		{
			focalT4IDs_patches_index = 0;
		}
		else
		{
			// Move focalT4IDs_patches_index to the next patch
			bool gotANewPatch = false;
			for (size_t i = focalT4IDs_patches_index + 1 ; i < focalT4IDs_patches.size() ; ++i)
			{
				if (focalT4IDs_patches[i] != focalT4IDs_patches[focalT4IDs_patches_index])
				{
					focalT4IDs_patches_index = i;
					gotANewPatch = true;
					break;
				}
			}
			assert(gotANewPatch);
		}	
		auto patch = focalT4IDs_patches[focalT4IDs_patches_index];


		// Compute stuff
		auto data = computeSmallSegmentsDiversityFromPaintedHaplotypes(allGatheredSegments[index_for_patch]);
		auto nbHaplotypesSampled = allGatheredSegments[index_for_patch].size();

		std::string generations_s = std::to_string(paintedGeneration == std::numeric_limits<int>::lowest() ? -GP->memory_nbGenerationsInBurnIn : paintedGeneration) + "\t" + std::to_string(observedGeneration) + "\t";
		//std::cout << "data.size() = " << data.size() << "\n";
		for (const auto& segment : data)
		{
			std::string s;
			if (SSP->SegDiversityFile_includeMainColor)
			{
				s = generations_s + std::to_string(segment.left) + "\t" + std::to_string(segment.right) + "\t" + std::to_string(patch) + "\t" + std::to_string(nbHaplotypesSampled) + "\t" + std::to_string(segment.nbColors) + "\t" + std::to_string(segment.colorHighestFrequency) + "\t" + std::to_string(segment.heterozygosity) + "\n";
			} else
			{
				s = generations_s + std::to_string(segment.left) + "\t" + std::to_string(segment.right) + "\t" + std::to_string(patch) + "\t" + std::to_string(nbHaplotypesSampled) + "\t" + std::to_string(segment.nbColors) + "\t" + std::to_string(segment.heterozygosity) + "\n";
			}
			file.write(s);
		}
	}

	file.close();


	/*
	std::pair<std::string, std::string> ss;
	{
		const auto data = computeSmallSegmentsFrequenciesFromPaintedHaplotypes(allGatheredSegments);

		ss.first.reserve(data.size() * (2 + 16 + 10));
		for (const auto& segment : data)
		{
			ss.first += "[" + std::to_string(segment.left) + " " + std::to_string(segment.right) + " " + std::to_string(segment.freq) + "]";
		}
	}

	{
		const auto data = computeLargeSegmentsFrequenciesFromPaintedHaplotypes(allGatheredSegments);
		
		ss.second.reserve(data.size() * (2 + 16 + 10));
		for (const auto& segment : data)
		{
			ss.second += "[" + std::to_string(segment.left) + " " + std::to_string(segment.right) + " " + std::to_string(segment.freq) + "]";
		}
	}
	*/
}

void T4TreeRec::writePaintedHaplotypes(const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches, const int paintedGeneration, const int observedGeneration, OutputFile& file) const
{
	//std::cout << "enters in writePaintedHaplotypes\n";
	assert(paintedGeneration < observedGeneration);
	assert(focalT4IDs_patches.size() == focalT4IDs.size());

	if (focalT4IDs.size() == 0) return;



	auto allGatheredSegments = computePaintedHaplotypes(paintedGeneration, focalT4IDs, focalT4IDs_patches);	




	int focalT4IDs_patches_index = -1;
	file.open();
	//std::cout << "allGatheredSegments.size() = " << allGatheredSegments.size() << "\n";
	for (size_t index_for_patch = 0 ; index_for_patch < allGatheredSegments.size() ; ++index_for_patch)
	{
		// Get patch 
		if (focalT4IDs_patches_index == -1)
		{
			focalT4IDs_patches_index = 0;
		}
		else
		{
			bool gotANewPatch = false;
			for (size_t i = focalT4IDs_patches_index + 1 ; i < focalT4IDs_patches.size() ; ++i)
			{
				if (focalT4IDs_patches[i] != focalT4IDs_patches[focalT4IDs_patches_index])
				{
					focalT4IDs_patches_index = i;
					gotANewPatch = true;
					break;
				}
			}
			assert(gotANewPatch);
		}	
		auto patch = focalT4IDs_patches[focalT4IDs_patches_index];

		// print stuff
		//std::cout << "allGatheredSegments[index_for_patch].size() = " << allGatheredSegments[index_for_patch].size() << "\n";
		for (size_t haplo_index = 0 ; haplo_index < allGatheredSegments[index_for_patch].size() ; ++haplo_index)
		{
			auto& gatheredSegments = allGatheredSegments[index_for_patch][haplo_index];


			// Write sampled haplotype info
			{
				auto I_index = haplo_index / 2;
				auto H_index = haplo_index % 2;
				std::string paintedGenerationString = paintedGeneration == std::numeric_limits<int>::lowest() ? std::to_string(-GP->memory_nbGenerationsInBurnIn) : std::to_string(paintedGeneration);
                std::string s(paintedGenerationString + "-" + std::to_string(observedGeneration) + " P" + std::to_string(patch) + " I" + std::to_string(I_index) + " H" + std::to_string(H_index) + ": ");
                //std::cout << "About to write " << s << "\n";
				file.write(s);
			}


			// Write segments
			assert(gatheredSegments[0].segment.left == 0);
			assert(gatheredSegments.back().segment.right == SSP->Gmap.T4_nbLoci);
			for (uint32_t i = 0 ; i < gatheredSegments.size() ; ++i)
			{
				auto& current  = gatheredSegments[i];
				if (i > 0)
				{
					auto& previous = gatheredSegments[i-1];
					assert(current.segment.left < current.segment.right);
					assert(previous.segment.left < previous.segment.right);
					assert(current.segment.left == previous.segment.right);
				}

				std::string s = "[" + std::to_string(current.segment.left) + " " + std::to_string(current.segment.right) + " P" + std::to_string(current.patch_index) + " I" + std::to_string(current.ind_ID) + "]";
				file.write(s);
			}

			// Write new line
			{
				std::string s = "\n";
				file.write(s);
			}
		}
	}
	file.close();
}


void T4TreeRec::computePaintedHaplotypes_exploreTree(HaplotypesContainer<HaplotypeOfSegments>& allSegments, const int paintedGeneration) const
{
	//this->print();
	for (int edge_index = edges.size()-1 ; edge_index >= 0 ; --edge_index )
	{
		auto& edge = edges[edge_index];

		assert(nodes.size() > edge.parent);
		if (nodes[edge.parent].generation >= paintedGeneration) 
		{
			// If the child is ancestor of the unique individual inserted in first place
			if (allSegments.doesAlreadyExist(edge.child))
			{
				// Get child
				HaplotypeOfSegments* childHaplotypeP = allSegments.getHaploP(edge.child);
				assert(childHaplotypeP != nullptr);

				// Get parent
				HaplotypeOfSegments* parentHaplotypeP;
				if (allSegments.doesAlreadyExist(edge.parent))
				{
					parentHaplotypeP = allSegments.getHaploP(edge.parent);
				} else
				{
					parentHaplotypeP = new HaplotypeOfSegments(nodes[edge.parent]);
					allSegments.insertHaploP(edge.parent, parentHaplotypeP);
				}

				/*
				if (nodes[edge.parent].generation == paintedGeneration)
				{
					std::cout << "Parent ID " << edge.parent << " existed at paintedGeneration = " << paintedGeneration << "\n";
				}
				*/

				/*
				std::cout << "-----------\nedge\n";
				edge.print();
				std::cout << "Before: parent ("<<edge.parent<<")\n";
				parentHaplotypeP->print();
				std::cout << "Before: child ("<<edge.child<<")\n";
				childHaplotypeP->print();
				*/
				
				

				//std::cout << edge_index << "\n";

				// Transfer segments (simple transfer means it does not "segmentize" and it does not update the tree)
				childHaplotypeP->sortAndMerge();
				parentHaplotypeP->simpleTransferSegmentsFromChild(*childHaplotypeP, edge, true);

				/*
				std::cout << "After: parent ("<<edge.parent<<")\n";
				parentHaplotypeP->print();
				std::cout << "After: child ("<<edge.child<<")\n";
				childHaplotypeP->print();
				*/
				


				// Delete if empty
				if (childHaplotypeP->size() == 0)
				{
					allSegments.deleteHaplo(edge.child);
				}
				if (parentHaplotypeP->size() == 0)
				{
					allSegments.deleteHaplo(edge.parent);
				}
			}
		}
	}
}


std::vector<std::vector<std::vector<T4TreeRec::PaintedSegment>>> T4TreeRec::computePaintedHaplotypes_gatherSegments(HaplotypesContainer<HaplotypeOfSegments>& allSegments, const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches) const
{
	std::vector<std::vector<std::vector<PaintedSegment>>> allPatchesGatheredSegments;

	assert(focalT4IDs_patches.size());
	assert(focalT4IDs_patches.size() == focalT4IDs.size());


	//////////////////////////////////////////////////////////////////
	// Allocate memory for allPatchesGatheredSegments and build map //
	//////////////////////////////////////////////////////////////////

	// allPatchesGatheredSegments[index_for_patch][index_for_haplo][segment_index]

	std::map<uint32_t, std::pair<size_t, size_t>> map;  //  allPatchesGatheredSegments[map[T4ID].first][map[T4ID].second]

	int index_for_patch = -1;
	int previousPatch = -1;
	int index1 = -1;
	int index2 = -1;
	for (size_t T4ID_absolute_index = 0 ; T4ID_absolute_index < focalT4IDs.size() ; ++T4ID_absolute_index)
	{
		assert(T4ID_absolute_index < focalT4IDs_patches.size());
		auto& patch = focalT4IDs_patches[T4ID_absolute_index];
		if (previousPatch != patch)
		{
			// new patch
			previousPatch = patch;
			++index_for_patch;
			allPatchesGatheredSegments.push_back({});
			index2 = 0;
			++index1;
		} else
		{
			++index2;
		}

		allPatchesGatheredSegments.back().push_back({});

		assert(index1 >= 0);
		assert(index2 >= 0);
		assert(allPatchesGatheredSegments.size() > index1);
		assert(allPatchesGatheredSegments[index1].size() > index2);

		assert(T4ID_absolute_index < focalT4IDs.size());
		auto focalT4ID = focalT4IDs[T4ID_absolute_index];
		map.insert( 
			std::pair<uint32_t, std::pair<size_t, size_t>>(
				focalT4ID,
				{
					index1,
					index2
				}
			)
		);
	}


	/////////////////////////////////////////////////////
	// Loop through allSegments to distribute segments //
	/////////////////////////////////////////////////////

	allSegments.iterator_restart();
	while (allSegments.iterator_isMore())
	{
		const auto info = allSegments.iterator_next();
		info.haploP->sortAndMerge();
		auto& segments = info.haploP->getSegmentsRef();

		for (auto& segment : segments)
		{
			//std::cout << "Testing segment.child is in the focals\n";
			assert(std::lower_bound(focalT4IDs.begin(),focalT4IDs.end(),segment.child) != focalT4IDs.end());

			/*
			std::cout << "Testing info.oldID lived in the generation 0\n";
			std::cout << "nodes["<<info.oldID<<"].generation = " << nodes[info.oldID].generation << "\n";
			std::cout << "nodes["<<info.oldID<<"] is associated here to segment = [" << segment.left << "-" << segment.right << "]" << std::endl;
			*/
			assert(nodes.size() > info.oldID);
			assert(nodes[info.oldID].generation == 0);

			
			auto& dest = map[segment.child];
			assert(allPatchesGatheredSegments.size() > dest.first);
			assert(allPatchesGatheredSegments[dest.first].size() > dest.second);
			allPatchesGatheredSegments[dest.first][dest.second].push_back({segment, info.haploP->getPatchIndex(), info.oldID}); // sort them later
		}
	}



	///////////////////////////////////////////////////////////////
	// Sort, merge and creat map of ID for the gathered segments //
	///////////////////////////////////////////////////////////////

	for (auto& patchGatheredSegments : allPatchesGatheredSegments)
	{
		for (auto& indGatheredSegments : patchGatheredSegments)
		{
			//// Sort the gathered segments 
			std::sort(
				indGatheredSegments.begin(),
				indGatheredSegments.end(), 
				[](const PaintedSegment & a, const PaintedSegment & b) -> bool
				{ 
				    return a.segment.left < b.segment.left; 
				}
			);

			//// Merge the gathered segments of same ancestors and create map of IDs
			if (indGatheredSegments.size() > 1)
			{	
				for (uint32_t i = 0 ; i < indGatheredSegments.size()-1 ; ++i)
				{
					auto l = i;
					while (indGatheredSegments.size() > l+1 && indGatheredSegments[l].ind_ID == indGatheredSegments[l+1].ind_ID && indGatheredSegments[l].segment.right == indGatheredSegments[l+1].segment.left)
					{
						++l;
					}
					if (l != i)
					{
						auto newRight = indGatheredSegments[l].segment.right;
						assert(newRight > indGatheredSegments[i].segment.right);
						indGatheredSegments[i].segment.right = newRight;
						indGatheredSegments.erase(indGatheredSegments.begin() + i + 1, indGatheredSegments.begin() + l + 1);
					}
				}
			}
		}
	}
		

	///////////////////////////
	// Free remaining memory //
	///////////////////////////

	allSegments.deleteAllHaplos();

	// Return
	return allPatchesGatheredSegments; // no newline at the end please!
}

std::vector<std::vector<std::vector<T4TreeRec::PaintedSegment>>> T4TreeRec::computePaintedHaplotypes(const int paintedGeneration, const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches) const
{
	#ifdef DEBUG
	std::cout << "Enters in std::string T4TreeRec::computePaintedHaplotype(int paintedGeneration) const\n\n\n\n";
	#endif
	HaplotypesContainer<HaplotypeOfSegments> allSegments;


	


	//this->print();
	/*
		It is a bottom-up approach. Starting with a unique individual and climbing up until it finds the ancestors
	*/


	//////////////////////////////
	// Insert unique individual //
	//////////////////////////////

	for (const auto& focalT4ID : focalT4IDs)
	{
		// insert unique individual
		assert(nodes[focalT4ID].generation == GP->CurrentGeneration);
		HaplotypeOfSegments* uniqueCurrentInd = new HaplotypeOfSegments({{0, SSP->Gmap.T4_nbLoci, (int)focalT4ID}}, nodes[focalT4ID]);
		uniqueCurrentInd->reset_alreadySegmentized(true);

		assert(focalT4ID < nodes.size());

		//std::cout << "focalT4ID = " << focalT4ID << "\n";

		allSegments.insertHaploP(
			focalT4ID,
			uniqueCurrentInd
		);
	}		


	//////////////////
	// Explore tree //
	//////////////////


	computePaintedHaplotypes_exploreTree(allSegments, paintedGeneration);


	/////////////////////////
	// Gather all segments //
	/////////////////////////

	//std::cout << "\n\n\n\nGathering segments...\n";

	return computePaintedHaplotypes_gatherSegments(allSegments, focalT4IDs, focalT4IDs_patches);
}


bool T4TreeRec::haveAllLociCoallesced_ifKnown() const
{
	assert(haveAllLociCoallesced_info == 'u' || haveAllLociCoallesced_info == 'y' || haveAllLociCoallesced_info == 'n');
	return haveAllLociCoallesced_info == 'y';
}


void T4TreeRec::simplify(Pop& pop)
{
	//std::cout << "Enter T4TreeRec::simplify!\n";
	//std::cout << "Entering simplify with " << nodes.size() << " nodes and " << edges.size() << " edges\n";
	#ifdef DEBUG
	//std::cout << "Entering simplify with " << nodes.size() << " nodes and " << edges.size() << " edges. Ratio of RAM usage = "<< (double)edges.size() * 4.0 / (double)nodes.size() <<"\n";
	std::cout << "Entering simplify with " << nodes.size() << " nodes and " << edges.size() << " edges\n";
	std::cout << "entering tree \n";
	this->print();
	#endif
	

	if (isAlreadyAsSimplifiedAsPossible) return;
	isAlreadyAsSimplifiedAsPossible = true;

	/*
	for (auto& elem : generationsToKeepInTheTree)
		std::cout << elem << " ";
	std::cout << "\n";
	*/
	lastGenerationSimplified = GP->CurrentGeneration;
	auto& Ei = edges;
	auto& Ni = nodes;
	NodeTable No;
	EdgeTable Eo;
	No.reserve(Ni.size() / 50);
	Eo.reserve(Ei.size() / 50);


	

	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Set up the segments of the current generation and the map of IDs for the current generation //
	/////////////////////////////////////////////////////////////////////////////////////////////////

	LastGenerationIDmap lastGenerationIDmap(Ni.size(), SSP->TotalpatchSize);

	HaplotypesContainer<HaplotypeOfSegments> allSegments;

	/*
	if (Ni.back().generation != GP->CurrentGeneration)
	{
		Ni.print();
		std::cout << "GP->CurrentGeneration = " << GP->CurrentGeneration << "\n";
		std::cout << "Ni.back().generation = " << Ni.back().generation << "\n";
	}
	*/

	assert(Ni.back().generation == GP->CurrentGeneration);
	//std::cout << "T4TreeRec::simplify: First RAM allocated\n";

	for ( int oldNodeID = Ni.size()-1 ; oldNodeID >= 0 && Ni[oldNodeID].generation == GP->CurrentGeneration ; --oldNodeID)
	{
		auto newNodeID = No.addNode(Ni[oldNodeID]);
		HaplotypeOfSegments* haploP = new HaplotypeOfSegments({Segment(0, SSP->Gmap.T4_nbLoci, newNodeID)}, Ni[oldNodeID]);
		haploP->reset_alreadySegmentized(true);
		haploP->set_newNodeID(newNodeID);
		allSegments.insertHaploP(
			oldNodeID,
			haploP
		);

		lastGenerationIDmap.add(oldNodeID, newNodeID);
	}
	//std::cout << "T4TreeRec::simplify: Current generation segments set\n";



	///////////////////////
	// Loop through tree //
	///////////////////////

	/*
	std::cout << "Ei.size() = " << Ei.size() << "\n";
	for (int edge_index = Ei.size() - 1 ; edge_index >= 0  ; --edge_index )
	{
		auto& edge = Ei[edge_index];
		std::cout << edge_index << ": {" << edge.left << " " << edge.right << " "<< edge.parent << " " << edge.child << "}\n";
	}*/

	
	for (int edge_index = Ei.size() - 1 ; edge_index >= 0  ; --edge_index )
	{
		auto& edge = Ei[edge_index];
		//std::cout << "--------------------\nedge: {" << edge.left << " " << edge.right << " "<< edge.parent << " " << edge.child << "}\n";

		if (allSegments.doesAlreadyExist(edge.child))
		{
			// Get child haplotype
			HaplotypeOfSegments* childHaplotypeP = allSegments.getHaploP(edge.child);
			assert(childHaplotypeP != nullptr);


			// Get parent haplotype
			HaplotypeOfSegments* parentHaplotypeP;

			if (allSegments.doesAlreadyExist(edge.parent))
			{
				parentHaplotypeP = allSegments.getHaploP(edge.parent);
			} else
			{
				parentHaplotypeP = new HaplotypeOfSegments(Ni[edge.parent]); // Note HaplotypeOfSegments does not know its own old ID but only its generation, its patch and its eventual future newID. Only allSegments knows the old ID
				allSegments.insertHaploP(edge.parent, parentHaplotypeP);
			}
			
			// Test if node must be kept in new tree.
			//std::cout << "ancestralGeneration = " << ancestralGeneration << "\n";
			//std::cout << "childHaplotypeP->getGeneration() = " << childHaplotypeP->getGeneration() << "\n";
			assert(childHaplotypeP->getGeneration() != ancestralGeneration);
			bool shouldNodeBeKeptIfItHasASegment = 
				generationsToKeepInTheTree.size()
				&&
				std::find(generationsToKeepInTheTree.begin(), generationsToKeepInTheTree.end(),childHaplotypeP->getGeneration()) != generationsToKeepInTheTree.end()
			;
			

			// do stuff

			//std::cout << "child oldID " << edge.child << " has segments:\n";

			/*
			{
				std::cout << "edge:"; edge.print(); 
				std::cout << "child before\n";
				std::cout << edge.child << " -> " << childHaplotypeP->getNewNodeID() << "\n";
				childHaplotypeP->print();
				std::cout << "------\n";
				std::cout << "parent before\n";
				std::cout << edge.parent << " -> " << parentHaplotypeP->getNewNodeID() << "\n";
				parentHaplotypeP->print();
				std::cout << "------\n";
			}
			*/

			parentHaplotypeP->transferSegmentsFromChild(*childHaplotypeP, edge, Eo, No, shouldNodeBeKeptIfItHasASegment);

			/*
			{
				std::cout << "child after\n";
				std::cout << edge.child << " -> " << childHaplotypeP->getNewNodeID() << "\n";
				childHaplotypeP->print();
				std::cout << "------\n";
				std::cout << "parent after\n";
				std::cout << edge.parent << " -> " << parentHaplotypeP->getNewNodeID() << "\n";
				parentHaplotypeP->print();
				std::cout << "------\n";
			}
			*/
			//std::cout << "newID is " << childHaplotypeP->getNewNodeID() << "\n"; 
			//std::cout << "----------------\n";

			// Remove child haplotype if empty (if it has no other parent)
			if (childHaplotypeP->size() == 0)
			{
				allSegments.deleteHaplo(edge.child);
			}
		} 
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Add ancestors to tree, reset ancestors, free remaining memory in allSegments and figure whether all loci have coalesced yet //
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	doStuffWithAncestors(allSegments, No, Eo);

	////////////////
	// Swap trees //
	////////////////
	No.swap(Ni);
	Eo.swap(Ei);


	//////////////////////
	// Reverse new tree //
	//////////////////////

	// It sets the newID into the so-called newnewID.
	Ni.reverse(); // That's lightning fast
	Ei.reverseAndMerge(Ni.size()); // That's slow


	assert(edges.back().child < nodes.size());
	//std::cout << "tree reversed\n";

	////////////////////////
	// Set new IDs to pop //
	////////////////////////

	//std::cout << "setting new IDs...\n";
	for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);
		for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			auto& ind = patch.getInd(ind_index);
			for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
			{
				auto& haplo = ind.getHaplo(haplo_index);
				auto newnewID = nodes.size() - lastGenerationIDmap.getNewID(haplo.T4ID) - 1;
				assert(newnewID >= 0 && newnewID < nodes.size());
				haplo.T4ID = newnewID;
			}
		}
	}

	#ifdef DEBUG
	std::cout << "Exiting simplify with " << nodes.size() << " nodes and " << edges.size() << " edges\nExiting tree:\n";
	this->print();
	std::cout << "###############################\n\n\n";
	#endif
	//std::cout << "Exiting simplify with " << nodes.size() << " nodes and " << edges.size() << " edges\nExiting tree:\n";
}




/*void T4TreeRec::drawSegments(std::vector<Segment>& segments)
{
	for (auto it = segments.cbegin(); it < segments.cend(); ++it)
	{
		auto l  = it->left;
		auto r = it->right;
		auto ID = it->newID;

		assert(l < r);

		std::cout << ID << ":";

		uint32_t i = 0;
		for ( ; i < l; ++i)
			std::cout << " ";
		for (; i < r; ++i)
			std::cout << "-";
		std::cout << "\n";
	}
}*/

void T4TreeRec::doStuffWithAncestors(HaplotypesContainer<HaplotypeOfSegments>& A, NodeTable& No, EdgeTable& Eo)
{
	///////////////////
	// Add ancestors //
	///////////////////
	#ifdef DEBUG
	std::cout << "tree just before doing stuff with ancestors\n";
	No.print();
	Eo.print();
	#endif

	uint32_t nbAncestorsAdded = 0;
	A.iterator_restart();
	while (A.iterator_isMore())
	{
		const auto info = A.iterator_next();
		auto& haploP = info.haploP;

		#ifdef DEBUG
		std::cout << "about to add ancestor with oldID " << info.oldID << "\n";
		if (haploP->getAlreadySegmentized())
		{
			std::cout << "Oops. ancestor newnewID undefined yet, newID undefined yet and oldID " << info.oldID << " has already been segmentized\n";
		}
		#endif
	
		assert(!haploP->getAlreadySegmentized());
		haploP->segmentize(No, Eo, true);
		++nbAncestorsAdded;
	}
	#ifdef DEBUG
	std::cout << "nbAncestorsAdded = " << nbAncestorsAdded << "\n";
	std::cout << "tree after adding ancestors\n";
	No.print();
	Eo.print();
	#endif

	////////////////////
	// Reset genetics //
	////////////////////
	std::vector<NodeGenetics> new_ancestorsGenetics;
	new_ancestorsGenetics.reserve(ancestorsGenetics.size());


	// To figure if all loci coalesced yet
	std::vector<Segment> allAncestorsSegment;

	uint32_t nbAncestors = 0; // just for assertion
	A.iterator_restart();
	while (A.iterator_isMore())
	{
		const auto info = A.iterator_next();
		auto& haploP = info.haploP;
		auto& oldID = info.oldID;
		auto& newID = info.newID;

		assert(nodes.size() > oldID);
		//std::cout << "nodes["<<oldID<<"].generation = " << nodes[oldID].generation << "\n";
		//std::cout << "ancestralGeneration = " << ancestralGeneration << "\n";
		assert(nodes[oldID].generation == ancestralGeneration);

		if (newID != -1)
		{
			assert(newID == haploP->getNewNodeID());
			auto newnewID = No.size() - newID - 1;
			assert(newnewID >= 0 && newnewID < No.size());

			++nbAncestors;

			auto oldAncestorIndex = std::lower_bound(
				ancestorsGenetics.begin(),
				ancestorsGenetics.end(),
				oldID
			) - ancestorsGenetics.begin();
			
			#ifdef DEBUG
			if (oldAncestorIndex >= ancestorsGenetics.size())
			{
				std::cout << "could not find oldID " << oldID << ". in ancestors. Ancestor IDs are:";
				for (auto& elem : ancestorsGenetics) std::cout << elem.getID() << " ";
				std::cout << "\n";
			}
			#endif
			assert(oldAncestorIndex < ancestorsGenetics.size());

			new_ancestorsGenetics.push_back(ancestorsGenetics[oldAncestorIndex]);
			new_ancestorsGenetics.back().resetID(newnewID);


			// To figure if all loci coalesced yet
			
			for (auto& segment : haploP->getSegmentsRef())
				allAncestorsSegment.push_back(segment);
			

			// free some memory but do not delete pointer yet so that I don't mess with the for loop.
			ancestorsGenetics[oldAncestorIndex].clear();
			ancestorsGenetics[oldAncestorIndex].shrink_to_fit();
		}			
	}	

	new_ancestorsGenetics.swap(ancestorsGenetics);
	std::sort(ancestorsGenetics.begin(), ancestorsGenetics.end());

	#ifdef DEBUG
	uint32_t nbAncestors2 = 0;
	for (uint32_t i = 0; i < No.size() ; ++i)
	{
		if (No[i].generation == ancestralGeneration) nbAncestors2++;
	}
	std::cout << "asserting number of ancestors: nbAncestors2 = " << nbAncestors2 << " ancestorsGenetics.size() = " << ancestorsGenetics.size() << " nbAncestors = " << nbAncestors << "\n";
	assert(nbAncestors2 == ancestorsGenetics.size());
	std::cout << "There are " << nbAncestors << " ancestors\n";
	#endif
	assert(nbAncestors == ancestorsGenetics.size());



	

	///////////////////////////
	// Free remaining memory //
	///////////////////////////
	A.deleteAllHaplos();


	//////////////////////////////////////
	// Figure if all loci coalesced yet //
	//////////////////////////////////////

	std::sort(allAncestorsSegment.begin(), allAncestorsSegment.end());

	/*
	std::cout << "printing all ancetors segments\n";
	for (auto& elem : allAncestorsSegment) elem.print();
	std::cout << "finished printing all ancetors segments\n";
	*/

	assert(allAncestorsSegment.front().left == 0);
	uint32_t to = 0;
	bool hasAllLociCoalesced = true;
	for (auto& ancestorSegments : allAncestorsSegment )
	{
		//ancestorSegments.print();
		assert(ancestorSegments.left <= to); // ensure no holes where loci would have somehow been lost
		if (ancestorSegments.left != to)
		{
			hasAllLociCoalesced = false;
			break;
		}
		to = ancestorSegments.right;
	}

	if (hasAllLociCoalesced)
	{
		assert(allAncestorsSegment.front().left == 0);
		assert(allAncestorsSegment.back().right == SSP->Gmap.T4_nbLoci);
		assert(to == SSP->Gmap.T4_nbLoci);
		haveAllLociCoallesced_info = 'y';
	} else
	{
		haveAllLociCoallesced_info = 'n';
	}

	/////////////////////////////////////////////////////////////////
	// Figure if one specific locus coalesced yet for killOnDemand //
	/////////////////////////////////////////////////////////////////
	if (killOnDemand_T4Locus != -1)
	{
		killOnDemand_isT4LocusFixed = true;
		bool hasOneSegmentHasAlreadyBeenFound = false;
		for (auto& ancestorSegments : allAncestorsSegment )
		{
			if (ancestorSegments.left > killOnDemand_T4Locus) break;
			if (killOnDemand_T4Locus < ancestorSegments.right && killOnDemand_T4Locus >= ancestorSegments.left)
			{
				if (hasOneSegmentHasAlreadyBeenFound)
				{
					killOnDemand_isT4LocusFixed = false;
					break;
				} else
				{
					hasOneSegmentHasAlreadyBeenFound = true;
				}
			}
			
		}

		killOnDemand_T4Locus = -1;
	}	
}


void T4TreeRec::setLocusForWhichFixationMustBeComputedAtTheNextSimplify(int locus)
{
	assert(killOnDemand_T4Locus == -1);
	killOnDemand_T4Locus = locus;
}

bool T4TreeRec::isLocusForWhichFixationHadToBeComputedFixed()
{
	return killOnDemand_isT4LocusFixed;
}

void T4TreeRec::print() const
{
	nodes.print();
	edges.print();
}


/*
template<typename INT> const std::vector<uint32_t>& getMutationsOfID(std::vector<std::vector<std::vector<uint32_t>>>& mutations, Pop& pop, INT ID)
{
	NodeGenetics node(ID);
	auto index = std::lower_bound(ancestorsGenetics.begin(), ancestorsGenetics.end(), node) - ancestorsGenetics.begin();
	assert(index >= 0 && index < ancestorsGenetics.size());
	return ancestorsGenetics[index].getMutations();
}
*/


template<typename INT> const std::vector<uint32_t>& T4TreeRec::getMutationsOfID(INT ID)
{
	NodeGenetics node(ID);
	auto index = std::lower_bound(ancestorsGenetics.begin(), ancestorsGenetics.end(), node) - ancestorsGenetics.begin();
	assert(index >= 0 && index < ancestorsGenetics.size());
	return ancestorsGenetics[index].getMutations();
}



void T4TreeRec::propagateMutationsForSegment(NodeGenetics& child, const NodeGenetics& parent, const uint32_t left, const uint32_t right)
{
	assert(left < right);
	//std::cout << "child.getGeneration() = " << child.getGeneration() << "\n";
	//std::cout << "parent.getGeneration() = " << parent.getGeneration() << "\n";

	assert(child.getGeneration() > parent.getGeneration());
	auto nbGenerationsInBetween = child.getGeneration() - parent.getGeneration();

	//whatHasBeenSetYet.push({left, right, -1});


	////////////////////////////////
	// Set segment as from parent //
	////////////////////////////////

	if (parent.getMutations().size())
	{
		// Get iterators
		auto insertFrom = std::lower_bound(parent.getMutations().cbegin(), parent.getMutations().cend(), left);    // lower only works
		auto insertTo = std::lower_bound(insertFrom, parent.getMutations().cend(), right);                   // lower only works
		auto startOfInsertion = std::lower_bound(child.getMutations().begin(), child.getMutations().end(), left); // lower or upper. Both should work

		// security -> Only one parent node can affect the a given segment
			
		if (startOfInsertion != child.getMutations().end()) 
		{
			
			// if something comes after the point of insertion then the next value must be greater or equal to right
			assert( *startOfInsertion >= left );

			if (startOfInsertion != (child.getMutations().end() - 1))
			{
				assert( *(startOfInsertion+1) >= right );
			}
		}

		// insert
		if (insertFrom != insertTo)
			child.getMutations().insert(startOfInsertion, insertFrom, insertTo);
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
        
        /*
        if (!(MutPosition>= left && MutPosition < right))
        {
        	std::cout << "right = " << right << "\n";
        	std::cout << "left = " << left << "\n";
        	std::cout << "MutPosition = " << MutPosition << "\n";
        }
        */
        
        //std::cout << MutPosition << " ";
        assert(MutPosition>= left && MutPosition < right);

        // Make the mutation
        //std::cout << "1309\n";
        child.mutateLocus(MutPosition);
        //std::cout << "1311\n";
    }
    //std::cout <<"\n";

    /*
    std::cout << "Childs genetics: ";
    for (auto& m : child.mutations) std::cout << m << " ";
    std::cout << "\n"; 
	*/
}


void T4TreeRec::deleteAll()
{
	edges.deleteAll();
	nodes.deleteAll();
	std::vector<NodeGenetics>().swap(ancestorsGenetics);
	std::vector<int>().swap(generationsToKeepInTheTree);
}


void T4TreeRec::PrintBinaryFile(OutputFile& file) const
{
	edges.PrintBinaryFile(file);
	nodes.PrintBinaryFile(file);

	file.writeBinary(ancestorsGenetics.size());
	for (auto& ancestorsGenetic : ancestorsGenetics)
	{
		ancestorsGenetic.PrintBinaryFile(file);
	}

	file.writeBinary(ancestralGeneration);
	file.writeBinary(lastGenerationSimplified);
	file.writeBinary(haveAllLociCoallesced_info);
	file.writeBinary(haveAllLociCoallesced_info);
}



void T4TreeRec::readFromBinaryFile(BinaryFileToRead& binfile)
{
	edges.readFromBinaryFile(binfile);
	nodes.readFromBinaryFile(binfile);

	{
		size_t ancestorsGeneticsSize;
		binfile.read(ancestorsGeneticsSize);
		ancestorsGenetics.resize(ancestorsGeneticsSize);

		for (auto& ancestorsGenetic : ancestorsGenetics)
		{
			ancestorsGenetic.readFromBinaryFile(binfile);
		}
	}

	binfile.read(ancestralGeneration);
	binfile.read(lastGenerationSimplified);
	binfile.read(haveAllLociCoallesced_info);
	binfile.read(haveAllLociCoallesced_info);
}
