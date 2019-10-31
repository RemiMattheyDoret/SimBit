template<typename vectorType, typename iteratorType>
ZipIterator<vectorType, iteratorType>::ZipIterator(const ZipIterator<vectorType, iteratorType>& other)
{
	this->haploP = other.haploP;
	this->flippedP = other.flippedP;
	this->haploEndP = other.haploEndP;
	this->flippedEndP = other.flippedEndP;
	this->skipPotentialValueFoundInBothHaploAndFlipped();
}


template<typename vectorType, typename iteratorType>
ZipIterator<vectorType, iteratorType>::ZipIterator(iteratorType a, iteratorType b, iteratorType c, iteratorType d)
{
	this->haploP = a;
	this->flippedP = b;
	this->haploEndP = c;
	this->flippedEndP = d;
	this->skipPotentialValueFoundInBothHaploAndFlipped();
}

template<typename vectorType, typename iteratorType>
ZipIterator<vectorType, iteratorType>::ZipIterator(){}

template<typename vectorType, typename iteratorType>
ZipIterator<vectorType, iteratorType>& ZipIterator<vectorType, iteratorType>::operator=(const ZipIterator<vectorType, iteratorType>&& other)
{
	this->haploP = other.haploP;
	this->flippedP = other.flippedP;
	this->haploEndP = other.haploEndP;
	this->flippedEndP = other.flippedEndP;
	this->skipPotentialValueFoundInBothHaploAndFlipped();
	return *this;
}

template<typename vectorType, typename iteratorType>
ZipIterator<vectorType, iteratorType>& ZipIterator<vectorType, iteratorType>::operator=(const ZipIterator<vectorType, iteratorType>& other)
{
	this->haploP = other.haploP;
	this->flippedP = other.flippedP;
	this->haploEndP = other.haploEndP;
	this->flippedEndP = other.flippedEndP;
	this->skipPotentialValueFoundInBothHaploAndFlipped();
	return *this;
}

template<typename vectorType, typename iteratorType>
void ZipIterator<vectorType, iteratorType>::skipPotentialValueFoundInBothHaploAndFlipped()
{
	while (flippedP != flippedEndP && haploP != haploEndP && *haploP == *flippedP)
	{
		++haploP;
		++flippedP;
	}	
}

template<typename vectorType, typename iteratorType>
unsigned int ZipIterator<vectorType, iteratorType>::operator*()
{

	if (flippedP != flippedEndP)
	{
		unsigned flippedPValue = *flippedP;
		if (haploP != haploEndP)
		{
			unsigned haploPValue = *haploP;
			if (haploPValue < flippedPValue)
			{
				return haploPValue;
			} else // if (*haploP > *flippedP) Should not have same value here thanks to 'skipPotentialValueFoundInBothHaploAndFlipped' who have been previously called
			{
				return flippedPValue;
			} 

		} else // if (haploP == haploEndP)
		{
			return flippedPValue;
		}
	} else // if (flippedP == flippedEndP)
	{
		return *haploP;
	}
}

template<typename vectorType, typename iteratorType>
bool operator==(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType>& rhs)
{
	return lhs.flippedP == rhs.flippedP && lhs.haploP == rhs.haploP;
}

template<typename vectorType, typename iteratorType>
bool operator!=(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType>& rhs)
{
	return lhs.haploP != rhs.haploP || lhs.flippedP != rhs.flippedP;
}

template<typename vectorType, typename iteratorType>
bool operator>(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType>& rhs)
{
	return lhs.flippedP > rhs.flippedP || lhs.haploP > rhs.haploP;
}

template<typename vectorType, typename iteratorType>
bool operator<(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType>& rhs)
{
	return lhs.flippedP < rhs.flippedP || lhs.haploP < rhs.haploP;
}


/*template<typename vectorType, typename iteratorType>
bool operator==(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs)
{
	return lhs.flippedP == rhs.flippedP && lhs.haploP == rhs.haploP;
}

template<typename vectorType, typename iteratorType>
bool operator!=(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs)
{
	return lhs.haploP != rhs.haploP || lhs.flippedP != rhs.flippedP;
}

template<typename vectorType, typename iteratorType>
bool operator>(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs)
{
	return lhs.flippedP > rhs.flippedP || lhs.haploP > rhs.haploP;
}

template<typename vectorType, typename iteratorType>
bool operator<(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs)
{
	return lhs.flippedP < rhs.flippedP || lhs.haploP < rhs.haploP;
}*/



template<typename vectorType, typename iteratorType>
ZipIterator<vectorType, iteratorType> ZipIterator<vectorType, iteratorType>::operator++()
{
	// std::cout << "flippedP = " << &(*flippedP) << "\n";
	// std::cout << "haploP = " << &(*haploP) << "\n";
	// std::cout << "flippedEndP = " << &(*flippedEndP) << "\n";
	// std::cout << "haploEndP = " << &(*haploEndP) << "\n";

	if (flippedP != flippedEndP)
	{
		unsigned flippedPValue = *flippedP;
		if (haploP != haploEndP)
		{
			unsigned haploPValue = *haploP;
			if (haploPValue < flippedPValue)
			{
				++haploP;
				this->skipPotentialValueFoundInBothHaploAndFlipped();
			} else // if (*haploP > *flippedP) Should not have same value here thanks to 'skipPotentialValueFoundInBothHaploAndFlipped' who have been previously called
			{
				++flippedP;
				this->skipPotentialValueFoundInBothHaploAndFlipped();
			} 

		} else // if (haploP == haploEndP)
		{
			++flippedP;
		}
	} else // if (flippedP == flippedEndP)
	{
		++haploP;
	}
	return *this;
}

/*template<typename vectorType, typename iteratorType>
void ZipIterator<vectorType, iteratorType>::skip()
{
	// std::cout << "flippedP = " << &(*flippedP) << "\n";
	// std::cout << "haploP = " << &(*haploP) << "\n";
	// std::cout << "flippedEndP = " << &(*flippedEndP) << "\n";
	// std::cout << "haploEndP = " << &(*haploEndP) << "\n";

	if (flippedP != flippedEndP)
	{
		if (haploP != haploEndP)
		{
			if (*haploP < *flippedP)
			{
				++haploP;
			} else // if (*haploP > *flippedP) Should not have same value here thanks to 'skipPotentialValueFoundInBothHaploAndFlipped' who have been previously called
			{
				++flippedP;
			} 
			this->skipPotentialValueFoundInBothHaploAndFlipped();

		} else // if (haploP == haploEndP)
		{
			++flippedP;
		}
	} else // if (flippedP == flippedEndP)
	{
			++haploP;
	}
}*/

/*ZipIterator<vectorType, iteratorType>& ZipIterator<vectorType, iteratorType>::operator++() // prefix
{
	std::cout << "call ZipIterator<vectorType, iteratorType>::operator++\n";
	if (haploP == haploEndP)
	{
		//std::cout << *flippedP << "\n";
		assert(flippedP != flippedEndP);
		++flippedP;
	} else if (flippedP == flippedEndP)
	{
		//std::cout << *haploP << "\n";
		assert(haploP != haploEndP);
		++haploP;
	} else
	{
		if (*haploP > *flippedP)
		{
			//std::cout << *flippedP << "\n";
			++flippedP;
		} else if (*haploP < *flippedP)
		{
			//std::cout << *haploP << "\n";
			++haploP;
		} else
		{
			do 
			{

				//std::cout << *flippedP << "\n";
				//std::cout << *haploP << "\n";
				++flippedP;
				++haploP;
			} while (flippedP != flippedEndP && haploP != haploEndP && *haploP == *flippedP);
		}
	}

	return *this;
}*/

template<typename vectorType, typename iteratorType>
bool ZipIterator<vectorType, iteratorType>::isMore()
{
	return haploP != haploEndP || flippedP != flippedEndP;
}

template<>
void ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator>::lower_bound(unsigned int target)
{
	haploP.lower_bound_FromThis(target);
	flippedP.lower_bound_FromThis(target);
	this->skipPotentialValueFoundInBothHaploAndFlipped();
}

template<>
void ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator>::lower_bound(unsigned int target)
{
	haploP   = std::lower_bound(haploP,   haploEndP,   target);
	flippedP = std::lower_bound(flippedP, flippedEndP, target);
	this->skipPotentialValueFoundInBothHaploAndFlipped();
}

template<typename vectorType, typename iteratorType>
void ZipIterator<vectorType, iteratorType>::upper_bound(unsigned int target)
{
	this->lower_bound(target);
	if (this->isMore()) // If end not reached
	{
		this->skip();
	}
}

template<typename vectorType, typename iteratorType>
iteratorType ZipIterator<vectorType, iteratorType>::get_haploP()
{
	return haploP;
}

template<typename vectorType, typename iteratorType>
iteratorType ZipIterator<vectorType, iteratorType>::get_flippedP()
{
	return flippedP;
}


