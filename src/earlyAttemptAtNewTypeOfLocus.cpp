
class Segment
{
	std::vector<bool> data;

public:
	double w;

	template<typename INT>
	bool operator[](INT i) {return data[i];}
	std::vector<bool>& getData() {return data;}

	Segment makeOffspring(std::vector<size_t>& newMuts);
};

class SegmentReference
{
	std::vector<bool> data;
	std::vector<Segment*> pointing;
};

Segment Segment::makeOffspring(std::vector<size_t>& newMuts)
{
	Problem because at each new mutation, we need to modify all segments that use this specific reference. I  need to use the logic that each segment use a vector<size_t> that flip the meaning of the reference.
}

class MapElement
{
	std::vector<SegmentReference> segmentReferences;
	std::vector<size_t> putaPoly;

	void removeNonPolymorphic();
};


void MapElement::removeNonPolymorphic()
{
	assert(segmentReferences.size());
	assert(segmentReferences.front().size());
	assert(segmentReferences.front().pointing.size());
	
	std::vector<bool> first = segmentReferences.front().pointing.front();
	std::vector<bool> stillPoly(false, first.size());
	assert(first.size() == putaPoly.size());


	//// Figure what is polymorphic
	for (auto& segmentReference : segmentReferences)
	{
		for (auto segmentPtr : segmentReference.pointing)
		{
			auto& segment = *segmentPtr;
			assert(segment.size() == first.size());
			bool areAllStillPoly = true;
			for (size_t i = 0 ; i < first.size() ; ++i)
			{
				if (!stillPoly[i])
				{
					if (segment[i] != first[i])
					{
						stillPoly[i] = true;
					}
				}

				if (!stillPoly[i])
				{
					areAllStillPoly = false;
				}
			}

			if (areAllStillPoly) return;
		}
	}


	//// Remove what is not polymorphic all segments
	for (auto& segmentReference : segmentReferences)
	{
		for (auto segmentPtr : segmentReference.pointing)
		{
			auto& segmentD = segmentPtr->getData();

			segmentD.erase(
				std::remove_if(
					segmentD.begin(),
					segmentD.end(),
					[&segmentD, &stillPoly](const auto& elem) { return  !stillPoly[&elem - &*segmentD.begin()]; }
				),
				segmentD.end()
			);
		}
	}

	//// Remove what is not polymorphic from putaPoly
	putaPoly.erase(
		std::remove_if(
			putaPoly.begin(),
			putaPoly.end(),
			[&putaPoly, &stillPoly](const auto& elem) { return  !stillPoly[&elem - &*putaPoly.begin()]; }
		),
		putaPoly.end()
	);
}
