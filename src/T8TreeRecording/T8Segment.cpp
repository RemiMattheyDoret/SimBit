


T8Segment* T8Segment::buildChild(T8Segment* parent, T8Segment* child)
{
	//std::cout << "building child " << child << " from parent " << parent << "\n";
	//assert(child != nullptr);
	//assert(parent != nullptr);
	assert(child != parent);

	child->parent = parent;
	child->nbChildren = 0;
	child->w = parent->w;
	//child->isSorted = true; // it contains no mutations for the moment so it can be considered sorted

	if (parent->nbChildren == std::numeric_limits<uint16_t>::max())
	{
		std::cout << "Oops, a single T8 segment has " << parent->nbChildren << " children! Tell Remi, he needs to use more than 16 bits to track the number of children. Could also be a bug, so please don't hesitate to report it to Remi.\n";
		abort();
	}

	++(parent->nbChildren);
	return child;
}

T8Segment* T8Segment::buildChild(std::vector<uint32_t>& genData, T8Segment* child)
{
	child->parent = nullptr;
	child->geneticData.swap(genData);
	genData.resize(0);
	if (SSP->T8_isSelection) child->forceComputeFitness(); else child->w = 1.0;
	child->nbChildren = 0;
	//child->isSorted = true; // the received genetic data must be sorted
	return child;
}

void T8Segment::forceComputeFitness()
{
	//assert(parent == nullptr);
	w = 1.0;
	if (SSP->T8_isSelection)
	{
		for (auto& mut : geneticData)
			w *= SSP->T8_FitnessEffects[mut];
		assert(w > 0.0);
	}
	//std::cout << "fitness computed = " << w << "\n";
}

template<typename INT> 
void T8Segment::mutate(const INT segmentIndex, bool mustBeSorted)
{
	//if (mustBeSorted) assert(std::is_sorted(geneticData.begin(), geneticData.end()));

	// sample number of mutations
	auto nbMuts = SSP->geneticSampler.get_T8_nbMuts(segmentIndex);
	if (nbMuts == 0) return;

	// reserve and resize
	auto newsize = geneticData.size() + nbMuts;
	if (newsize > geneticData.capacity())
		geneticData.reserve(newsize);
	
	else if (geneticData.capacity() > 20*newsize)
	{
		// Shrink if it needs to
		geneticData = {geneticData.begin(), geneticData.end()};
		geneticData.reserve(newsize);
	}


	// for each mutation
	for (size_t mut = 0 ; mut < nbMuts ; ++mut )
	{
		// sample mutation position
		auto position = SSP->geneticSampler.get_T8_mutationPosition(segmentIndex);

		// Make mutation
		if (mustBeSorted)
		{
			auto it = std::lower_bound(geneticData.begin(), geneticData.end(), position);
			geneticData.insert(it, position);
		} else
		{
			geneticData.push_back(position);
		}
			
		// Compute fitness - Can eventually use parental fitness
		if (SSP->T8_isSelection) w *= SSP->T8_FitnessEffects[position];
	}

	//if (mustBeSorted) assert(std::is_sorted(geneticData.begin(), geneticData.end()));
}

void T8Segment::clear_shrink()
{
	std::vector<uint32_t>().swap(geneticData);
	parent = nullptr;
}

void T8Segment::shrink_to_fit()
{
	geneticData.shrink_to_fit();
}

void T8Segment::clear()
{
	geneticData.resize(0);
	parent = nullptr;
}


