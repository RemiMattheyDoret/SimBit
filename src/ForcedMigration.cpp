
void ForcedMigration::readInput(InputReader& input)
{
	assert(SSP != nullptr);

	if (input.PeakNextElementString() == "default")
	{
		input.skipElement();
		return;
	}

	while (input.IsThereMoreToRead())
	{
		auto type = input.GetNextElementString();
		if (type != "move" && type != "copy" && type != "m" && type != "c") 
		{
			std::cout << "For option --forcedMigration, expected either keyword 'move' (or 'm') or 'copy' (or 'c') but received '" << type << "' instead.\n";
			abort();
		}
		if ((type == "move" || type == "m") && SSP->fecundityForFitnessOfOne == -1)
		{
			std::cout << "For option --forcedMigration, you asked to 'move' individuals. Individuals can only be moved if the fecundity differs from -1 (that is if the patch sizes can differ from the patch carrying capacity).\n";
			abort();
		}

		auto nbInds = input.GetNextElementInt();

		if (nbInds < 0)
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying a negative number of individuals. Received '" << type << " " << nbInds <<  "'.\n";
			abort();
		}

		/*
		auto fromtostring = input.GetNextElementString();
		if (fromtostring != "fromto" && fromtostring != "fromTo" && fromtostring != "FromTo" && fromtostring != "ft" && fromtostring != "FT")
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying " << nbInds << ". Expected keyword 'fromto' (or equivalent such as 'ft' or 'fromTo' for examples) but received '" << fromtostring << " instead.\n";
			abort();
		}
		*/

		auto patch_from = input.GetNextElementInt();
		auto patch_to = input.GetNextElementInt();

		auto atstring = input.GetNextElementString();
		if (atstring != "at")
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying " << nbInds << " from patch " << patch_from << " to patch " << patch_to << ". Expected keyword 'at' but received '" << atstring << " instead. Note that the patch numbers have not been checked at that moment yet.\n";
			abort();
		}

		auto isEarlyString = input.GetNextElementString();
		if (isEarlyString != "early" && isEarlyString != "late" && isEarlyString != "e" && isEarlyString != "l")
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying " << nbInds << " from patch " << patch_from << " to patch " << patch_to << ". Expected keyword 'early' or 'late' after keyword 'at' but received '" << isEarlyString << " instead. Note that the patch numbers have not been checked at that moment yet.\n";
			abort();
		}

		auto generation = input.GetNextElementInt();

		if (generation < 0 || generation > GP->nbGenerations)
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying " << nbInds << " from patch " << patch_from << " to patch " << patch_to << " at generation "<< generation << ". The generation received is either negative or greater than the total number of generations (which you set to "<<GP->nbGenerations<<").\n";
			abort();
		}

		int generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), generation) - GP->__GenerationChange.begin() - 1;

		assert(GP->__PatchNumber.size() > generation_index);

		if (patch_from < 0 || patch_to < 0)
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying " << nbInds << " from patch " << patch_from << " to patch " << patch_to << " at generation "<< generation << ". Either the 'patch_from' or the 'patch_to' is negative.\n";
			abort();
		}

		if (GP->__PatchNumber[generation_index] <= patch_from || GP->__PatchNumber[generation_index] <= patch_to)
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying " << nbInds << " from patch " << patch_from << " to patch " << patch_to << " at generation "<< generation << ". At generation " << generation << " there are only " << GP->__PatchNumber[generation_index] << " patches. The patch indicated is zero based conuting and must therefore be lower than the total number of patches.\n";
			abort();
		}

		if (patch_from == patch_to)
		{
			std::cout << "For option --forcedMigration, you asked for moving or copying " << nbInds << " from patch " << patch_from << " to patch " << patch_to << " at generation "<< generation << ". Source and destination patch must differ. Cannot move or copy individuals from a patch to the same patch.\n";
			abort();
		}

		if (isEarlyString == "early" || isEarlyString == "e")
		{
			earlytimes.push_back(generation);
			from_atEarlyTimes.push_back(patch_from);
			to_atEarlyTimes.push_back(patch_to);
			nbInds_atEarlyTimes.push_back(nbInds);
			isCopy_atEarlyTimes.push_back(type == "copy" || type == "c");
		} else
		{
			latetimes.push_back(generation);
			from_atLateTimes.push_back(patch_from);
			to_atLateTimes.push_back(patch_to);
			nbInds_atLateTimes.push_back(nbInds);
			isCopy_atLateTimes.push_back(type == "copy" || type == "c");
		}		
	}

	if (earlytimes.size())
	{
		auto idx = reverse_stable_sort_indexes(earlytimes);
		reorder(earlytimes, idx);
		reorder(from_atEarlyTimes, idx);
		reorder(to_atEarlyTimes, idx);
		reorder(nbInds_atEarlyTimes, idx);
		reorder(isCopy_atEarlyTimes, idx);
	}

	if (latetimes.size())
	{
		auto idx = reverse_stable_sort_indexes(latetimes);
		reorder(latetimes, idx);
		reorder(from_atLateTimes, idx);
		reorder(to_atLateTimes, idx);
		reorder(nbInds_atLateTimes, idx);
		reorder(isCopy_atLateTimes, idx);
	}

	assert(earlytimes.size() == from_atEarlyTimes.size());
	assert(earlytimes.size() == to_atEarlyTimes.size());
	assert(earlytimes.size() == nbInds_atEarlyTimes.size());
	assert(earlytimes.size() == isCopy_atEarlyTimes.size());

	assert(latetimes.size() == from_atLateTimes.size());
	assert(latetimes.size() == to_atLateTimes.size());
	assert(latetimes.size() == nbInds_atLateTimes.size());
	assert(latetimes.size() == isCopy_atLateTimes.size());
}



void ForcedMigration::forceMigrationIfNeeded(Pop& pop, bool isEarly)
{
	std::deque<Individual> migrantPool_inds;
	std::deque<size_t> migrantPool_dest;

	auto& times = isEarly ? earlytimes : latetimes ;

	while (times.size() && times.back() == GP->CurrentGeneration)
	{
		auto& from_atTimes = isEarly ? from_atEarlyTimes : from_atLateTimes;
		auto& to_atTimes = isEarly ? to_atEarlyTimes : to_atLateTimes;
		auto& nbInds_atTimes = isEarly ? nbInds_atEarlyTimes : nbInds_atLateTimes;
		auto& isCopy_atTimes = isEarly ? isCopy_atEarlyTimes : isCopy_atLateTimes;


		assert(SSP != nullptr);
		times.pop_back();
		auto from = from_atTimes.back(); from_atTimes.pop_back();
		auto to = to_atTimes.back(); to_atTimes.pop_back();
		auto nbInds = nbInds_atTimes.back(); nbInds_atTimes.pop_back();
		auto isCopy = isCopy_atTimes.back(); isCopy_atTimes.pop_back();

		assert(from < GP->PatchNumber);
		assert(to   < GP->PatchNumber);

		auto& patchFrom = pop.getPatch(from);
		auto& patchTo   = pop.getPatch(to);

		
		nbInds = nbInds < SSP->patchCapacity[to] ? nbInds : SSP->patchCapacity[to] ;
		if (nbInds == 0) continue;

		std::vector<size_t> ind_indices_from;

		// Create ind_indices_from by sampling without replacement
		{
			if (nbInds == 1)
			{
				ind_indices_from.push_back(GP->rngw.uniform_int_distribution(SSP->patchSize[from]));
			} else
			{
				assert(nbInds > 1);
				std::vector<size_t> x;
				x.reserve(SSP->patchSize[from]);				

				ind_indices_from.reserve(nbInds);				
				for (size_t i = 0 ; i < nbInds ; ++i)
				{
					if (x.size() == 0)
					{
						for (size_t i = 0 ; i < SSP->patchSize[from] ; ++i)
							x.push_back(i);
					}

					size_t j;
					if (x.size() == 1)					
					{
						j = 0;
						ind_indices_from.push_back(x[j]);
					}
					else
					{
						j = GP->rngw.uniform_int_distribution(x.size());
						ind_indices_from.push_back(x[j]);
						std::swap(x[j], x.back());
					}
					x.pop_back();
				}
			}
			assert(ind_indices_from.size() == nbInds);
		}
				
		//std::cout << "nbInds = " << nbInds << "\n";
		size_t ind_index_to = SSP->patchSize[to];
		for (size_t fake_ind_index_from = 0 ; fake_ind_index_from < nbInds ; ++fake_ind_index_from)
		{
			if (SSP->patchSize[from] && SSP->patchCapacity[to])
			{
				auto ind_index_from = ind_indices_from[fake_ind_index_from];
				//std::cout << "Sampled " << ind_index_from << " from " << SSP->patchSize[from] << "\n";
				//assert(ind_index_from >= 0 && ind_index_from < SSP->patchSize[from]);

				if (ind_index_to == SSP->patchCapacity[to])
					ind_index_to = 0;

				if (!isCopy)
				{
					patchTo.getInd(ind_index_to).swap(patchFrom.getInd(ind_index_from));
				} else
				{
					patchTo.getInd(ind_index_to) = patchFrom.getInd(ind_index_from);
				}
				//std::cout << "ind_index_to = " << ind_index_to << "\n";
				//std::cout << "SSP->patchCapacity["<<to<<"] = " << SSP->patchCapacity[to] << "\n";
				//std::cout << "SSP->patchSize["<<to<<"] = " << SSP->patchSize[to] << "\n";
				if (SSP->patchSize[to] < SSP->patchCapacity[to])
				{
					assert(ind_index_to == SSP->patchSize[to]);
					++(SSP->patchSize[to]);
				}
				++ind_index_to;
					
				
				if (!isCopy)
				{
					assert(SSP->fecundityForFitnessOfOne != -1.0);
					if (ind_index_from != SSP->patchSize[from]-1)
						patchFrom.getInd(ind_index_from).swap(patchFrom.getInd(SSP->patchSize[from]-1));

					--(SSP->patchSize[from]);
				}
			}
		}
	}
}
