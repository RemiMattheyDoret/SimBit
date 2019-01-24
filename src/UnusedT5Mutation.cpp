T5Mutation::T5Mutation(int l)
:locus(l), generation(GP->currentGeneration), ID(currentID)
{
	currentID++;
}


T5Mutation::operator==(T5Mutation other)
{
	if (this->ID != other.ID || this->generation != other.generation || this->locus != other.locus)
	{
		return false;
	}
	return true;
}


T5Mutation::operator==(int otherLocus)
{
	return this->locus == otherLocus;
}
