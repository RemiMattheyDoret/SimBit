T8TreeRecording::T8TreeRecording()
{
	//nbGenerationsSinceLastSimplify = 0;
}

void T8TreeRecording::printT8IDs()
{
	std::cout << "Parent pop has T8IDs up to (non-included) " << lastGenerationSegments.size() << "\n";
	std::cout << "Offspring pop has T8IDs up to (non-included) " << currentGenerationSegments.size() << "\n";
}

void T8TreeRecording::printT8IDs(Pop& pop)
{
	std::cout << "print T8 IDs: ";
	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);
		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			auto& ind = patch.getInd(ind_index);
			auto& haplo0 = ind.getHaplo(0);
			auto& haplo1 = ind.getHaplo(1);

			std::cout << haplo0.T8ID << " " << haplo1.T8ID << " ";
		}
	}
	std::cout << "\n";
}

void T8TreeRecording::initialize(Pop& pop)
{
	/*
	size_t mapIndex = 0;
	hash.resize(SSP->Gmap.T8_nbLoci);
	for (size_t locus = 0 ; locus < SSP->Gmap.T8_nbLoci; ++locus)
	{
		if (locus > SSP->T8_map[mapIndex]) ++mapIndex;
		//std::cout << "["<<locus << "-" << mapIndex << "]\n";
		hash[locus] = mapIndex;
	}
	*/


	// Assuming al loci starts fix for 0 (wildtype) variant

	assert(GP->startAtGeneration == GP->CurrentGeneration);
	// LastTimeMutationsWerePlaced
	lastTimeMutationsWerePlaced = GP->startAtGeneration;

	// For shirnk to fit
	lastGenerationMemPoolsWereShrunkToFit = GP->CurrentGeneration;

	// Initialize memory pool
	size_t totalTotalPatchCapacity=0;
	for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; speciesIndex++)
    {
        totalTotalPatchCapacity += allParameters.SSPs[speciesIndex].TotalpatchCapacity;
    }

    assert(SSP->T8_map.size());
    memPool_segments.setMultiplier(1.1);
	if (memPool_segments.size() == 0) memPool_segments.resize(SSP->T8_map.size()*4*totalTotalPatchCapacity);


	// Create first nodes and assign them to the right haplotypes
	assert(SSP->TotalpatchSize > 0);
	assert(GP->PatchNumber > 0);
	assert(SSP->patchSize.size() == GP->PatchNumber );

	currentGenerationSegments.resize(2 * SSP->TotalpatchSize);
	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);

		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			auto& ind = patch.getInd(ind_index);

			// Get haplotypes
			auto& haplo0 = ind.getHaplo(0);
			auto& haplo1 = ind.getHaplo(1);

			// set IDs
			haplo0.T8ID = currentGenerationIDcounter++;
			haplo1.T8ID = currentGenerationIDcounter++;

			//std::cout << "initialize ID at " << haplo0.T8ID << "\n";
			//std::cout << "initialize ID at " << haplo1.T8ID << "\n";

			assert(haplo0.T8ID < 2 * SSP->TotalpatchSize);
			assert(haplo1.T8ID < 2 * SSP->TotalpatchSize);


			// set fitnesses
			haplo0.setW_T8(1.0);
			haplo1.setW_T8(1.0);

			// Create all segments
			assert(SSP->T8_map.size());
			for (size_t segmentIndex = 0 ; segmentIndex < SSP->T8_map.size() ; ++segmentIndex)
			{
				auto node0 = memPool_segments.NEW(true); // for haplo0
				auto node1 = memPool_segments.NEW(true); // for haplo1
				assert(node0->nbChildren == 0);
				assert(node1->nbChildren == 0);
				assert(node0->parent == nullptr);
				assert(node1->parent == nullptr);
				assert(currentGenerationSegments.size() > haplo0.T8ID);
				assert(currentGenerationSegments.size() > haplo1.T8ID);
				currentGenerationSegments[haplo0.T8ID].push_back(node0);
				currentGenerationSegments[haplo1.T8ID].push_back(node1);
			}
		}
	}
}

void T8TreeRecording::deleteTree()
{
	//std::cout << "Just before deleting T8 tree there is " << totalNbNodes << " nodes.\n";
	prune(currentGenerationSegments);
	lastGenerationSegments.clear(); lastGenerationSegments.shrink_to_fit();
	currentGenerationSegments.clear(); currentGenerationSegments.shrink_to_fit();
}

T8TreeRecording::~T8TreeRecording()
{
	deleteTree();
}

void T8TreeRecording::deleteNode(T8Segment* node)
{
	if (node->parent == nullptr)
	{
		memPool_segments.DELETE(node, false);
	} else
	{
		memPool_segments.DELETE(node, true);
	}
}

void T8TreeRecording::prune(std::vector<std::vector<T8Segment*>>& childrenSegments)
{
	//std::cout << "Begin pruning...";
	// Prune does not simplify, in the sense that it does not remove nodes that serves no purpose. propagateMutations does that.

	size_t ID = 0;
	for (auto& segments : childrenSegments)
	{
		//std::cout << "pruning ID " << ID << "\n";
		++ID;
		for (auto segment : segments)
		{
			auto focal = segment;
			assert(focal != nullptr);
			while (focal->nbChildren == 0)
			{
				// Get parent
				auto parent = focal->parent;

				// delete focal
				deleteNode(focal);

				// Break if focal was an ancestor
				if (parent == nullptr) break;
		
				// Remove one child from parent
				--(parent->nbChildren);

				// set focal to parent
				focal = parent;

			}
		}
	}
	//std::cout << "End pruning\n";
}

void T8TreeRecording::announceEndOfGeneration()
{
	//std::cout << "## Before pruning\n";
	//memPool_segments.computeSomeStats();
	//if (GP->CurrentGeneration % 500 == 0) memPool_segments.computeSomeStats();
	prune(lastGenerationSegments);
	/*
	std::cout << "## After pruning\n";
	memPool_segments.computeSomeStats();
	std::cout << "-------\n";
	*/
	
	/*
	++nbGenerationsSinceLastSimplify;
	if (nbGenerationsSinceLastSimplify == SSP->T8_simplifyEveryNGenerations)
	{
		simplify();
		nbGenerationsSinceLastSimplify = 0;
	}*/
}

void T8TreeRecording::announceStartOfGeneration()
{
	currentGenerationIDcounter = 0;
	lastGenerationSegments.swap(currentGenerationSegments);

	if (currentGenerationSegments.size() != 2*SSP->TotalpatchSize) 
		currentGenerationSegments.resize(2*SSP->TotalpatchSize);

	if (lastGenerationMemPoolsWereShrunkToFit < GP->CurrentGeneration)
	{
		if (lastGenerationMemPoolsWereShrunkToFit == GP->startAtGeneration)
		{
			if (
				(GP->CurrentGeneration - lastGenerationMemPoolsWereShrunkToFit) 
				% 
				((int)std::max(10, (int)(SSP->TotalpatchSize/10)))
				== 
				0
			)
			{
				memPool_segments.moveBigToSmall();
				memPool_segments.emptySmallMemory();
				memPool_segments.shrink_all_to_fit();
				lastGenerationMemPoolsWereShrunkToFit = GP->CurrentGeneration;
			}
		} else
		{
			if (
				(GP->CurrentGeneration - lastGenerationMemPoolsWereShrunkToFit) 
				% 
				((int)std::max(500, (int)(SSP->TotalpatchSize/2.2)))
				== 
				0
			)
			{
				memPool_segments.moveBigToSmall();
				memPool_segments.emptySmallMemory();
				memPool_segments.shrink_all_to_fit();
				lastGenerationMemPoolsWereShrunkToFit = GP->CurrentGeneration;
			}
		}
	}
}


std::pair<T8ID_type, fitnesstype> T8TreeRecording::addOffspring(T8ID_type motherID, T8ID_type fatherID, std::vector<uint32_t>& recPositions)
{
	/*
		The terms father and mother are confusing because the parents are haplotypes. They obviously have no gender.
	*/

	// Deal with IDs
	/*
	if (!(lastGenerationSegments.size() > fatherID && lastGenerationSegments.size() > motherID))
	{
		std::cout << "fatherID = " << fatherID << "\n";
		std::cout << "motherID = " << motherID << "\n";
		std::cout << "lastGenerationSegments.size() = " << lastGenerationSegments.size() << "\n";
	}
	std::cout << "fatherID = " << fatherID << "\n";
	std::cout << "motherID = " << motherID << "\n";
	*/
	
	/*
	std::cout << "Enters addOffspring\n";
	std::cout << "fatherID = " << fatherID << "\n";
	std::cout << "motherID = " << motherID << "\n";
	std::cout << "recPositions = [";
	for (auto& e : recPositions) std::cout << e << " ";
	std::cout << "]\n";

	std::cout << "SSP->T8_map = [";
	for (auto& e : SSP->T8_map) std::cout << e << " ";
	std::cout << "]\n";
	*/

	//recPositions = {recPositions.back()};

	// Test IDs
	//assert(lastGenerationSegments.size() > motherID);

	// Take second parent segments if needed
	/*if (recPositions.size() > 1)
	{
		assert(lastGenerationSegments.size() > fatherID);
	}*/

	//std::cout << "SSP->T8_map: ";
	//printVector(SSP->T8_map);
		

	auto childID = currentGenerationIDcounter++;
	//assert(currentGenerationSegments.size() > childID);
	auto& childSegments  = currentGenerationSegments[childID];
	if (childSegments.size() != SSP->T8_map.size()) childSegments.resize(SSP->T8_map.size());
	

	// Loop through each segment and copy appropriately
	fitnesstype haplotypeFitness = 1.0;
	bool isMother = true;
	size_t recIndex = 0;
	size_t beginningOfSegment = 0;

	//std::cout << "SSP->T8_map: "; printVector(SSP->T8_map);
	//std::cout << "recPositions: "; printVector(recPositions);
	
	std::vector<bool> didRecombine; didRecombine.reserve(SSP->T8_map.size());
	
	for (size_t segmentIndex = 0 ; segmentIndex < SSP->T8_map.size() ; ++segmentIndex)
	{
		//std::cout << "segmentIndex = " << segmentIndex << " | recIndex = " << recIndex << " | beginningOfSegment = " <<beginningOfSegment<<"\n";
		size_t endOfSegment = SSP->T8_map[segmentIndex]; // end is included.
		assert(recIndex < recPositions.size()); // recPositions always end with INT_MAX

		// Ignore recombination in between blocks
		if (recPositions[recIndex] == endOfSegment)
		{
			// recombination in between segments
			isMother = !isMother;
			++recIndex;
			assert(recPositions.size() > recIndex);
		}

		// Check if recombination happend within this segment or not
		if (motherID != fatherID && recPositions[recIndex] < endOfSegment)
		{
			didRecombine.push_back(true);
			//// Recombination happened within this segment
			//std::cout << "Recombination at position "<<recPositions[recIndex]<<" happened within the segment "<< beginningOfSegment << "-"<< endOfSegment <<"\n";

			// Propagate mutation to both parent haplotypes
			auto& mother = lastGenerationSegments[motherID][segmentIndex];
			auto& father = lastGenerationSegments[fatherID][segmentIndex];
			
			
			propagateMutations(mother);
			propagateMutations(father);
			

			// Create recombined geneticData
			auto from = beginningOfSegment;
			auto to = recPositions[recIndex];


			
			
			// data have already been sorted
			while ( true )
			{
				std::vector<uint32_t>::iterator itFrom;
				std::vector<uint32_t>::iterator itTo;
				if (isMother)
				{
					itFrom = std::lower_bound(mother->geneticData.begin(), mother->geneticData.end(), from);
					itTo = std::lower_bound(itFrom, mother->geneticData.end(), to);
					//motherIt = itTo;
				} else
				{
					itFrom = std::lower_bound(father->geneticData.begin(), father->geneticData.end(), from);
					itTo = std::lower_bound(itFrom, father->geneticData.end(), to);
					//fatherIt = itTo;
				}
				assert(itFrom <= itTo);
				if (itTo > itFrom)
					geneticDataAtRecPlaces.insert(geneticDataAtRecPlaces.end(), itFrom, itTo);

				if (to == endOfSegment) break;

				++recIndex;
				from = to;
				to = recPositions[recIndex] < endOfSegment ? recPositions[recIndex] : endOfSegment ;

				isMother = !isMother;
			}
			
			

			// Get isMother and recIndex ready for next segment

			// Create the child segment
			childSegments[segmentIndex] = T8Segment::buildChild(
				geneticDataAtRecPlaces,  // swaps it and resize the new 'geneticDataAtRecPlaces', hence no need to geneticDataAtRecPlaces.resize(0)
				memPool_segments.NEW(false)
			);

			// Mutate this new segment
			childSegments[segmentIndex]->mutate(segmentIndex, true);

		} else
		{
			didRecombine.push_back(false);
			//// segment transfered without recombination
			childSegments[segmentIndex] = T8Segment::buildChild(
				isMother ? lastGenerationSegments[motherID][segmentIndex] : lastGenerationSegments[fatherID][segmentIndex],
				memPool_segments.NEW(true)
			);

			// Mutate this new segment
			childSegments[segmentIndex]->mutate(segmentIndex, false);
		}


		// compute fitness
		haplotypeFitness *= childSegments[segmentIndex]->w;

		// iterate
		beginningOfSegment = endOfSegment;
	}

/*
	auto dist = std::poisson_distribution<size_t>(SSP->T8_Total_Mutation_rate);
	auto nbMuts = dist(GP->rngw.getRNG());
	//std::cout << "nbMuts = " << nbMuts << "\n";
	for (size_t mut = 0 ; mut < nbMuts ; ++mut)
	{
		auto position = GP->rngw.uniform_int_distribution(SSP->Gmap.T8_nbLoci);
		assert(position < SSP->Gmap.T8_nbLoci);
		assert(hash[position] < childSegments.size());
		//std::cout << "hash["<<position<<"] = " << hash[position] << "\n";
		auto node = childSegments[hash[position]];
		if (didRecombine[hash[position]])
		{
			node->geneticData.insert(std::lower_bound(node->geneticData.begin(), node->geneticData.end(), position), position);
		} else
		{
			node->geneticData.push_back(position);
		}
			
		node->w *= SSP->T8_FitnessEffects[position];
	}

	haplotypeFitness = 1.0;
	for (size_t segmentIndex = 0 ; segmentIndex < SSP->T8_map.size() ; ++segmentIndex)
		haplotypeFitness *= childSegments[segmentIndex]->w;
*/

	//assert(haplotypeFitness >= 0.0);

	//std::cout << "haplotypeFitness = " << haplotypeFitness << "\n";
	//std::cout << "Exits addOffspring\n";
	return {childID, haplotypeFitness};
}


void T8TreeRecording::moveGeneticDataFromTo(T8Segment* giver, T8Segment* receiver)
{
	if (giver->geneticData.size() > receiver->geneticData.size())
	{
		giver->geneticData.insert(
			giver->geneticData.end(),
			std::make_move_iterator(receiver->geneticData.begin()),
			std::make_move_iterator(receiver->geneticData.end())
		);

		giver->geneticData.swap(receiver->geneticData);
		
	} else
	{
		receiver->geneticData.insert(
			receiver->geneticData.end(),
			std::make_move_iterator(giver->geneticData.begin()),
			std::make_move_iterator(giver->geneticData.end())
		);
	}
}


void T8TreeRecording::copyGeneticDataFromTo(T8Segment* giver, T8Segment* receiver)
{
	receiver->geneticData.insert(
		receiver->geneticData.end(),
		giver->geneticData.begin(),
		giver->geneticData.end()
	);
}

void T8TreeRecording::moveSortedGeneticDataFromTo(T8Segment* giver, T8Segment* receiver)
{
	size_t middle;
	if (giver->geneticData.size() > receiver->geneticData.size())
	{
		middle = giver->geneticData.size();
		giver->geneticData.insert(
			giver->geneticData.end(),
			receiver->geneticData.begin(),
			receiver->geneticData.end()
		);
		
		giver->geneticData.swap(receiver->geneticData);
	} else
	{
		
		middle = receiver->geneticData.size();
		receiver->geneticData.insert(
			receiver->geneticData.end(),
			giver->geneticData.begin(),
			giver->geneticData.end()
		);
		
	}

	std::inplace_merge(
		receiver->geneticData.begin(),
		receiver->geneticData.begin() + middle,
		receiver->geneticData.end()
	);
}


void T8TreeRecording::copySortedGeneticDataFromTo(T8Segment* giver, T8Segment* receiver)
{

	/*	
	auto middle = receiver->geneticData.size();
	receiver->geneticData.insert(
		receiver->geneticData.end(),
		giver->geneticData.begin(),
		giver->geneticData.end()
	);
		

	std::inplace_merge(
		receiver->geneticData.begin(),
		receiver->geneticData.begin() + middle,
		receiver->geneticData.end()
	);
	*/
	
	staticGeneticData.resize(0);
		
	std::merge(
		giver->geneticData.begin(), 
		giver->geneticData.end(),
		receiver->geneticData.begin(),
		receiver->geneticData.end(),
		std::back_inserter(staticGeneticData)
	);

	staticGeneticData.swap(receiver->geneticData);
}



void T8TreeRecording::propagateMutations(T8Segment* focal)
{
	//std::cout << "propagating mutations at generation " << GP->CurrentGeneration << "\n";
	// This function can probably be implemented in a smarter way. Maybe with a K-way merge algorithm or by sorting at the end instead of inplace merge.

	// add all new mutations along the lines and remove those that lose parenthood. Also simplify the lineage (remove nodes that contribute to only a single child)
	
	assert(focal != nullptr);

	//std::cout << "focal = " << focal << "\n";

	auto child = focal;
	auto parent = child->parent;
	//std::cout << "Propagating mutations for child " << child << " with parent " << parent << "\n";
	if (parent == nullptr) return;

	// Move focal to other memoryPool
	memPool_segments.makeBig(focal);


	/*
		Two solutions:
			
			1) I could walk the way up to build the current individual (and remove individuals with a single child from the tree). When doing that I must be careful not to make memory leaks by failing to delete nodes.
			
			2) I could walk up the tree (and eventually remove individuals with a single child from the tree) and build the vector of pointers to child. Then, I could accumulate mutations down to each nodes. This technique is slower but it will reduce the tree, hence making future call faster. But each node will have more genetic stuff so it might end up taking more RAM. Not sure.

	*/

	assert(SSP->T8_propagationMethod == 1 || SSP->T8_propagationMethod == 2);
	

	bool hasPassedFirstCAyet = false;
	while (parent->parent != nullptr)
	{
		assert(parent->nbChildren);


		if (parent->nbChildren == 1)
		{
			T8Segment* grandParent = parent->parent;

			if (hasPassedFirstCAyet)
			{
				assert(focal != child);
				copyGeneticDataFromTo(parent, focal);
			} else
			{
				assert(focal == child);
			}
			moveGeneticDataFromTo(parent, child);

			parent->nbChildren = 0;
			deleteNode(parent);

			child->parent = grandParent;
			parent = grandParent;

		} else
		{

			if (hasPassedFirstCAyet)
			{
				assert(focal != child);
				copyGeneticDataFromTo(parent, focal);

				child = parent;
				parent = parent->parent;
			} else
			{
				assert(focal == child);
				copyGeneticDataFromTo(parent, child);

				--(parent->nbChildren);

				child = parent;
				parent = parent->parent;

				hasPassedFirstCAyet = true;
			}
		}
	}


	// Deal with ancestor
	std::sort(focal->geneticData.begin(), focal->geneticData.end());
	assert(parent->parent == nullptr);
	if (parent->nbChildren == 1)
	{
		T8Segment* grandParent = parent->parent;

		if (hasPassedFirstCAyet)
		{
			assert(focal != child);
			copySortedGeneticDataFromTo(parent, focal);
			std::sort(child->geneticData.begin(), child->geneticData.end());
			moveSortedGeneticDataFromTo(parent, child);
		} else
		{
			assert(focal == child);
			moveSortedGeneticDataFromTo(parent, focal);
		}
			
		parent->nbChildren = 0;
		deleteNode(parent);

		child->parent = grandParent;
		//parent = grandParent;

	} else
	{
		if (hasPassedFirstCAyet)
		{
			assert(focal != child);
			copySortedGeneticDataFromTo(parent, focal);

			//child = parent;
			//parent = parent->parent;
		} else
		{
			assert(focal == child);
			copySortedGeneticDataFromTo(parent, child);

			--(parent->nbChildren);

			//child = parent;
			//parent = parent->parent;

			//hasPassedFirstCAyet = true;
		}
	}
		

	focal->parent = nullptr;
}


void T8TreeRecording::placeMutations()
{
	if (lastTimeMutationsWerePlaced == GP->CurrentGeneration) return;
	lastTimeMutationsWerePlaced = GP->CurrentGeneration;
	for (auto& segments : currentGenerationSegments)
	{
		assert(segments.size() == SSP->T8_map.size());
		for (size_t segmentIndex = 0 ; segmentIndex < segments.size() ; ++segmentIndex)
		{
			auto segment = segments[segmentIndex];
			propagateMutations(segment);

			auto& v = segment->geneticData; // for convenience

			// Remove duplicate so that 1 means any number of mutations. That's an annoying reduction in information for the output but I guess I'll be happy with that for the moment.
			if (SSP->T8_WhenToSortData == '0') // if not sorted during or at the end of propagation
			{
				std::sort(v.begin(), v.end());
				v.erase( std::unique( v.begin(), v.end() ), v.end() ); // That's questionable
			}
				

			// Assertions
			{
				if (v.size())
				{
					if (segmentIndex == 0)
					{
						assert(v.front() >= 0);	
					} else
					{
						assert(v.front() >= SSP->T8_map[segmentIndex-1]);
					}
					assert(v.back() < SSP->T8_map[segmentIndex]);
				}
			}
		}
	}
}

std::vector<std::vector<std::vector<uint32_t>>> T8TreeRecording::getMutationsData(Pop& pop, bool shouldDeleteTree)
{
	placeMutations();
	
	std::vector<std::vector<std::vector<uint32_t>>> r;
	r.resize(GP->PatchNumber);
	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		auto& patch = pop.getPatch(patch_index);
		r[patch_index].resize(SSP->patchSize[patch_index] * 2);
		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			auto& ind = patch.getInd(ind_index);

			// Haplotype 0
			{
				auto& segments = currentGenerationSegments[ind.getHaplo(0).T8ID];
				assert(segments.size() == SSP->T8_map.size());
				assert(segments.size());
				if (segments.size() == 1)
				{
					if (shouldDeleteTree)
						r[patch_index][ind_index * 2].swap(segments.back()->geneticData);
					else
						r[patch_index][ind_index * 2] = segments.back()->geneticData;
				} else
				{
					for (auto& segment : segments)
					{
						r[patch_index][ind_index * 2].insert(
							r[patch_index][ind_index * 2].end(),
							segment->geneticData.begin(),
							segment->geneticData.end()
						);

						if (shouldDeleteTree)
							std::vector<uint32_t>().swap(segment->geneticData);
					}
				}
			}

			// Haplotype 1
			{
				auto& segments = currentGenerationSegments[ind.getHaplo(1).T8ID];
				assert(segments.size() == SSP->T8_map.size());
				assert(segments.size());
				if (segments.size() == 1)
				{
					if (shouldDeleteTree)
						r[patch_index][ind_index * 2 + 1].swap(segments.back()->geneticData);
					else
						r[patch_index][ind_index * 2 + 1] = segments.back()->geneticData;
				} else
				{
					for (auto& segment : segments)
					{
						r[patch_index][ind_index * 2 + 1].insert(
							r[patch_index][ind_index * 2 + 1].end(),
							segment->geneticData.begin(),
							segment->geneticData.end()
						);

						if (shouldDeleteTree)
							std::vector<uint32_t>().swap(segment->geneticData);
					}
				}
			}
		}	
	}

	if (shouldDeleteTree)
	{
		deleteTree();
	}

	/*
	std::cout << "Data:\n\t";
	for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
	{
		for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
		{
			for (auto& mut : r[patch_index][ind_index * 2])
				std::cout << mut << " ";
			std::cout << "\n\t";

			for (auto& mut : r[patch_index][ind_index * 2 + 1])
				std::cout << mut << " ";
			std::cout << "\n\t";
		}
	}
	*/
	
	return r;
}


void T8TreeRecording::computeFrequenciesFromRawData(std::vector<std::vector<std::vector<uint32_t>>>& data, std::vector<std::map<uint32_t, double>>& freqs, bool shouldClearData)
{
	if (freqs.size() == 0) // Otherwise assume frequencies are computed and are fine
	{			
		// Compute frequencies
		assert(data.size() == GP->PatchNumber);
		freqs.resize(GP->PatchNumber);
		for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
		{
			assert(data[patch_index].size() == 2*SSP->patchSize[patch_index]);

			for (size_t haplo_index = 0 ; haplo_index < 2*SSP->patchSize[patch_index] ; ++haplo_index)
			{
				for ( auto& mut : data[patch_index][haplo_index] )
				{
					++freqs[patch_index][mut];
				}
			}
		}
	} else
	{
		assert(data.size() == 0);
	}

	// get rid of data
	if (shouldClearData)
		std::vector<std::vector<std::vector<uint32_t>>>().swap(data);
}
