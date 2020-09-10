
void Genealogy::printBufferIfIsTime()
{
	if (isItTime)
	{
		assert(outputWriter.isFile(genealogy));
    auto& file = outputWriter.get_OutputFiles(genealogy)[0];
		buffer += "\n";
		file.open();
		file.write(buffer);
		file.close();
		buffer = "";
		//std::cout << "buffer printed to " << file.getPath() << "\n";
	}
}

void Genealogy::startNewGeneration()
{
	assert(buffer.size() == 0);
	assert(isItTime);
	buffer = "G_" + std::to_string(GP->CurrentGeneration);
}

void Genealogy::addOffspringIfIsTime(int offPatch, int offIndex, int motherPatch, int motherIndex, int fatherPatch, int fatherIndex)
{
	if (isItTime)
	{
		buffer += whitespace + P + std::to_string(offPatch) + I + std::to_string(offIndex) + underscore;
		buffer += P + std::to_string(motherPatch) + I + std::to_string(motherIndex) + underscore;
		buffer += P + std::to_string(fatherPatch) + I + std::to_string(fatherIndex);
		//std::cout << "buffer = "<< buffer << "\n";
	}
}


void Genealogy::setupBothParentTables(std::vector<std::vector<bool>>& parentsForNext, std::vector<std::vector<bool>>& parents)
{
	assert(SSP->maxEverpatchCapacity.size() == GP->maxEverPatchNumber);
	for (int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; patch_index++)
	{
		std::vector<bool> v(SSP->maxEverpatchCapacity[patch_index], false);
		parents.push_back(v);
		parentsForNext.push_back(v);
		
		//Assertions
		/*
		assert(parents[patch_index].size() == SSP->maxEverpatchCapacity[patch_index]);
		assert(parentsForNext[patch_index].size() == SSP->maxEverpatchCapacity[patch_index]);
		for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
		{
			assert(!parents[patch_index][ind_index]);
			assert(!parentsForNext[patch_index][ind_index]);
		}
		*/
	}
}

void Genealogy::resetParentsToFalse(std::vector<std::vector<bool>>& Parents)
{
	assert(Parents.size() == GP->maxEverPatchNumber);
	for (int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; patch_index++)
	{
		assert(Parents[patch_index].size() == SSP->maxEverpatchCapacity[patch_index]);
		for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
		{
			Parents[patch_index][ind_index] = false;
		}
	}
}

void Genealogy::getFields(std::vector<std::string>& fields, std::string& token, std::string filePath)
{
	//std::cout << token << "\n";
	fields.reserve(3);
	tokenize(fields, token, '_');

  	// Assertions
  	if (fields.size() != 3)
  	{
  		std::cout << "In Genealogy::getFields. The token "<<token<<" of the file " << filePath << " contained " << fields.size() - 1 << " underscores instead of the two expected. This is either an internal problem or the user has modified the temporary files.\n";
  		abort();
  	}
  	if (fields[1].at(0) != 'P')
  	{
  		std::cout << "In Genealogy::getFields. A parent ID does not start with in the file " << filePath << " This is either an internal problem or the user has modified the temporary files.\n";
  		abort();
  	}
  	if (fields[2].at(0) != 'P')
  	{
  		std::cout << "In Genealogy::getFields. A parent ID does not start with in the file " << filePath << " This is either an internal problem or the user has modified the temporary files.\n";
  		abort();
  	}
}

void Genealogy::setParents(std::vector<std::string>& fields, std::vector<std::vector<bool>>& parents)
{
		//std::cout << "Begin of setParents parents[0][1] = " << parents[0][1] << "\n";

		// first parent
  	{
  		int fieldIndex = 1;
  		int Ipos = fields[fieldIndex].find('I');
    	assert(Ipos < fields.size());
    	//std::cout << "fields[fieldIndex] = " << fields[fieldIndex] << " : " << fields[fieldIndex].substr(1,Ipos-1) << " : " << fields[fieldIndex].substr(Ipos+1) << "\n";
    	int patch_index = std::stoi(fields[fieldIndex].substr(1,Ipos-1));
    	int ind_index   = std::stoi(fields[fieldIndex].substr(Ipos+1));

    	/*if (parents[patch_index][ind_index])
    	{
				std::cout << "P" << patch_index << "I" << ind_index << " (aka. "<<fields[fieldIndex]<<") was already parent\n";    		
    	}*/

      parents[patch_index][ind_index] = true;

      //std::cout << "P" << patch_index << "I" << ind_index << " (aka. "<<fields[fieldIndex]<<") will be a parent\n";
  	}

      // second parent
      {
  		int fieldIndex = 2;
  		int Ipos = fields[fieldIndex].find('I');
    	assert(Ipos < fields.size());
    	//std::cout << "fields[fieldIndex] = " << fields[fieldIndex] << " : " << fields[fieldIndex].substr(1,Ipos-1) << " : " << fields[fieldIndex].substr(Ipos+1) << "\n";
    	int patch_index = std::stoi(fields[fieldIndex].substr(1,Ipos-1));
    	int ind_index   = std::stoi(fields[fieldIndex].substr(Ipos+1));

    	/*if (parents[patch_index][ind_index])
    	{
				std::cout << "P" << patch_index << "I" << ind_index << " (aka. "<<fields[fieldIndex]<<") was already parent\n";    		
    	}*/

      parents[patch_index][ind_index] = true;
      //std::cout << "P" << patch_index << "I" << ind_index << " (aka. "<<fields[fieldIndex]<<") will be a parent\n";
  	}

  	//std::cout << "End of setParents parents[0][1] = " << parents[0][1] << "\n";
}


void Genealogy::getParentsIDsFromFile(OutputFile& file, std::vector<std::vector<bool>>& parents)
{
	std::string content;
	file.openAndReadLine(content, GP->CurrentGeneration);

	std::stringstream lineStream(content);
  std::string token;
  lineStream >> token; // Remove generation token
  while(lineStream >> token) // lineStream >> token
  {
  	// get the three fields
		std::vector<std::string> fields;
  	getFields(fields, token, file.getPath());

  	setParents(fields, parents);
  }
}

bool Genealogy::isIndividualParent(std::string& description, std::vector<std::vector<bool>>& parents)
{
	int Ipos = description.find("I");
  assert(Ipos < description.size());
  //std::cout << "description = " << description << " : " << description.substr(1,Ipos-1) << " : " << description.substr(Ipos+1) << "\n";
  int patch_index = std::stoi(description.substr(1,Ipos-1));
  int ind_index   = std::stoi(description.substr(Ipos+1));
  /*if (parents[patch_index][ind_index])
  {
  	std::cout << description << " is a parent\n";
  } else
  {
  	std::cout << description << " is not a parent\n";
  }*/
	return parents[patch_index][ind_index];
}

void Genealogy::rewriteFile(OutputFile& file, int generation, std::vector<std::vector<bool>>& parents, std::vector<std::vector<bool>>& parentsForNext)
{
	//std::cout << "Generation " << generation << "\n\n";

	std::string content;
	file.openAndReadLine(content, generation);

	std::stringstream lineStream(content);
  std::string token;
  std::string newLine;
  lineStream >> newLine; // Get Generation

  //int NBINDS = 0;
  //int NBPARENTS = 0;

  while(lineStream >> token)
  {
  	// get the three fields
  	//std::cout << "token = " << token << "\n";
		std::vector<std::string> fields;
  	getFields(fields, token, file.getPath(generation));

  	//NBINDS++;
  	if (isIndividualParent(fields[0], parents))
  	{
  		//NBPARENTS++;
  		newLine += whitespace + token;
  		setParents(fields, parentsForNext);
  	}
  }
  newLine += "\n";

  //std::cout << "rewriting file " << file.getPath(generation) << " | NBINDS = " << NBINDS << " | NBPARENTS = " << NBPARENTS << "\n";

  file.clearContentAndLeaveOpen(generation);
  file.write(newLine);
  file.close();	
}

bool Genealogy::isTimeToCoalesce()
{
	assert(coalesceGenealogyFrequency > 0);
	if (GP->CurrentGeneration > generationFrom)
	{
		if (GP->CurrentGeneration == generationTo)
		{
			return true;
		}
		if(
			(GP->CurrentGeneration < generationTo)
			&&
			((GP->CurrentGeneration - generationFrom) % coalesceGenealogyFrequency) == 0
		)
		{
			return true;
		}
	}
		
	return false;
}

void Genealogy::coalesce(OutputFile& file)
{
	assert(GP->CurrentGeneration <= generationTo);

	std::vector<std::vector<bool>>  parentsForNext;
	std::vector<std::vector<bool>>  parents;

	setupBothParentTables(parents, parentsForNext);

	for (int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; patch_index++)
	{
		assert(parents[patch_index].size() == SSP->maxEverpatchCapacity[patch_index]);
		assert(parentsForNext[patch_index].size() == SSP->maxEverpatchCapacity[patch_index]);
		for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
		{
			assert(!parents[patch_index][ind_index]);
			assert(!parentsForNext[patch_index][ind_index]);
		}
	}

	getParentsIDsFromFile(file, parents);

	//writeOutWhoIsParent(parents);

	for (int generation = GP->CurrentGeneration - 1 ; generation >= generationFrom ; generation--)
	{
		rewriteFile(file, generation, parents, parentsForNext);
		
		parents = parentsForNext;
		if (generation != generationFrom)
		{
			resetParentsToFalse(parentsForNext);
		}
  	//writeOutWhoIsParent(parents);
	}
}


void Genealogy::writeOutWhoIsParent(std::vector<std::vector<bool>>& parents)
{
	for (int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; patch_index++)
	{
		assert(parents[patch_index].size() == SSP->maxEverpatchCapacity[patch_index]);
		for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
		{
			if (parents[patch_index][ind_index])
			{
				std::cout << "confirm that P" << patch_index << "I" << ind_index << " is parent\n";
			} else
			{
				std::cout << "confirm that P" << patch_index << "I" << ind_index << " is not parent\n";
			}
		}
	}
}

bool Genealogy::isTimeToMerge()
{
	return GP->CurrentGeneration == generationTo;
}

void Genealogy::removeFilesAtStart(OutputFile& file)
{
	for (int i = 0 ; i < (generationTo - generationFrom + 1) ;i++)
	{
		int generation = generationFrom + i;
		remove(file.getPath(generation).c_str());
	}
}

void Genealogy::mergeFiles(OutputFile& file)
{
	std::vector<std::string> listOfFiles(generationTo - generationFrom + 1);
	for (int i = 0 ; i < listOfFiles.size() ;i++)
	{
		int generation = generationFrom + i;
		listOfFiles[i] = file.getPath(generation);
	}
	file.mergeFiles(listOfFiles);

	for (std::string& p : listOfFiles)
		remove(p.c_str());
}


bool Genealogy::isTime()
{
	//bool b = GP->CurrentGeneration >= generationFrom && GP->CurrentGeneration <= generationTo;
	//std::cout << "willItEverBeTime = " << willItEverBeTime << " GP->CurrentGeneration = " << GP->CurrentGeneration << " generationFrom = " << generationFrom << " generationTo = " << generationTo << " second test = " << b << "\n";
	if (!willItEverBeTime)
	{
		isItTime = false;
		return isItTime;
	}
	if (GP->CurrentGeneration >= generationFrom && GP->CurrentGeneration <= generationTo)
	{
		isItTime = true;
		//std::cout << "set isItTime to true\n";
		return isItTime;
	}

	isItTime = false;
	return isItTime;
}

bool Genealogy::isCoalesce()
{
	return coalesceGenealogyFrequency > 0;
}



void Genealogy::setGenealogyToNothing()
{
	generationFrom = -1;
	generationTo = -1;
	isItTime = false;
	willItEverBeTime = false;
	coalesceGenealogyFrequency = -1;
}

void Genealogy::setGenealogyTimes(std::vector<int> times)
{
	if (times.size() != 2)
	{
		std::cout << "Error from Genealogy::setGenealogyTimes. In input for --genealogy_file, it appears that you have given more than two time points. Only two time points can be given for this option, the generation from which SimBit starts to compute the genealogy and the generation at which SimBit stops computing the genealogy.\n";
		abort();
	}
	if (times[0] > times[1])
	{
		std::cout << "Error in Genealogy::setGenealogyTimes. the first time received is bigger than the second time received. This is an internal error but you might want to check your input for --genealogy_file anyway!\n";
		abort();
	}
	generationFrom = times[0];
	generationTo = times[1];
	willItEverBeTime = true;
}

void Genealogy::setcoalesceGenealogyFrequency(int x)
{
	coalesceGenealogyFrequency = x;
}


