

std::vector<T4TreeRec::Segment>::iterator T4TreeRec::HaplotypeOfSegments::begin(){return segments.begin();}
std::vector<T4TreeRec::Segment>::iterator T4TreeRec::HaplotypeOfSegments::end(){return segments.end();}

void T4TreeRec::HaplotypeOfSegments::assertOrderingAndNoOverlap()
{
	for (uint32_t i = 0 ; i < segments.size() ; ++i)
	{
		auto& segment = segments[i];
		assert(segment.left < segment.right);
		if (i)
		{
			assert(segment.left >= segments[i].right);
		}
	}
}

void T4TreeRec::HaplotypeOfSegments::removeMovingSequence(T4TreeRec::Segment& movingSegment, uint32_t& childSegmentIndex)
{
	assert(alreadySegmentized);
	assert(childSegmentIndex < segments.size());
	auto& segment = segments[childSegmentIndex];
	assert(movingSegment.left >= segment.left && movingSegment.right <= segment.right);
	
	if (movingSegment.left > segment.left && movingSegment.right < segment.right)
	{
		/*
			childSegment    ------------
			MSeq               -----

			Need to duplicate
		*/

		// First insert new childSegment
		// Watch out this invalidates the reference Segment
		segments.insert(
			segments.begin() + childSegmentIndex,
			{
				segment.left,
				movingSegment.left,
				segment.child
			}		
		);

		// Second correct other childSegment
		segments[childSegmentIndex+1].left = movingSegment.right;

		++childSegmentIndex;

	} else
	{
		/*
			childSegment    -------
			MSeq            -------

			childSegment    -------
			MSeq            ----

			childSegment    -------
			MSeq               ----

			No need to duplicate
		*/

		if (segment.left < movingSegment.left)
		{
			/*
				childSegment    -------
				MSeq               ----
			*/
			#ifdef DEBUG
			assert(movingSegment.right == segment.right);
			#endif
			segment.right = movingSegment.left;
		} else if (segment.right > movingSegment.right)
		{
			/*
				childSegment       -------
				MSeq               ----
			*/
			#ifdef DEBUG
			assert(movingSegment.left == segment.left);
			#endif
			segment.left = movingSegment.right;
		} else
		{
			/*
				childSegment       -------
				MSeq               -------
			*/
			#ifdef DEBUG
			assert(movingSegment.left == segment.left);
			assert(movingSegment.right == segment.right);
			#endif

			segments.erase(segments.begin() + childSegmentIndex); // This invalidates the reference segments
			--childSegmentIndex;
		}
	}
}

void T4TreeRec::HaplotypeOfSegments::print() const
{
	
    for (auto& elem : this->segments)
    {
        std::cout << "{" << elem.left << ", " << elem.right << " (" << elem.child << ")} ";
    }
    std::cout << "\n";
    
    /*
    for (auto& elem : this->segments)
    {
        for (uint32_t i = 0 ; i < elem.left ; ++i)
            std::cout << " ";
        for (uint32_t i = 0 ; i < elem.right - elem.left ; ++i)
            std::cout << "-";
        std::cout << "  " << elem.child << " {" << elem.left << " " << elem.right << "}\n";           
    }
    std::cout << "\n\n\n";  
    */
}


T4TreeRec::Segment T4TreeRec::HaplotypeOfSegments::getMovingSegment(const Edge& edge, T4TreeRec::Segment& childSegment)
{
	auto& MSleft = std::max(edge.left, childSegment.left);
	auto& MSright = std::min(edge.right, childSegment.right);

	#ifdef DEBUG
	assert(MSleft < MSright);
	#endif

	return T4TreeRec::Segment(
		MSleft,
		MSright,
		childSegment.child
	);
}

void T4TreeRec::HaplotypeOfSegments::sortAndMerge()
{
	alreadySegmentized = true;

	std::sort(segments.begin(), segments.end());
	
	for (uint32_t i = 0 ; i < segments.size()-1 ; ++i)
	{
		auto l = i;
		while (segments.size() < l+1 && segments[l].child == segments[l+1].child && segments[l].right == segments[l+1].left)
		{
			++l;
		}
	
		if (l != i)
		{
			auto newRight = segments[l].right;
			assert(newRight > segments[i].right);
			segments[i].right = newRight;
			#ifdef DEBUG
			std::cout << "i = " << i << "\n";
			std::cout << "l = " << l << "\n";
			#endif
	
			segments.erase(segments.begin() + i + 1, segments.begin() + l + 1);
		}
	}
}

void T4TreeRec::HaplotypeOfSegments::segmentize(T4TreeRec::NodeTable& nodeTable, T4TreeRec::EdgeTable& edgeTable, bool shouldNodeBeKeptIfItHasASegment)
{
	assert(!alreadySegmentized);
	alreadySegmentized = true;

	#ifdef DEBUG
	if (shouldNodeBeKeptIfItHasASegment)
	{
		std::cout << "Node that should be kept:\n";
		for (auto const & seg : segments)
			seg.print();
		std::cout << "End of segments of Node that should be kept\n";
	}
	#endif
	


	// If no segments
	if (segments.size() == 0) // This scenario is possible as haplotype pointers are created before knowing whether any segment will be passed through a specific edge
	{
		return;
	}




	// If one segment
	else if (segments.size() == 1)
	{
		if (shouldNodeBeKeptIfItHasASegment)
		{
			// Special case, where the segment must go through current individual
			if (this->newNodeID == -1)
			{
				this->newNodeID = nodeTable.addNode(myNode);
			} 
			assert(this->newNodeID > segments[0].child);
			edgeTable.addEdge({segments[0].left, segments[0].right, this->newNodeID, segments[0].child});
			segments = {{segments[0].left, segments[0].right, this->newNodeID}};
		}

		return;
	}


	/*
		-------
		   ---------------
		          ---
		               ------
								---
								  ---
								        ---  
		At any time, the next event is always either
			- previous right
			- next left  (recursive)
			- next right (recursive)
			- my right

	*/


	// generate all transition points
    std::vector<T4TreeRec::PointForFindingOverlap> points;
    for (auto const & seg : segments)
    {
        points.push_back({seg.left, true, seg.child});
        points.push_back({seg.right, false, seg.child});
    }

    // sort transition points
    std::sort(points.begin(), points.end(), 
      [](const T4TreeRec::PointForFindingOverlap& a, const T4TreeRec::PointForFindingOverlap& b) { return a.location < b.location; });

    std::vector<T4TreeRec::Segment> B;
    B.reserve(segments.size());

    // initialize overlaps
    std::vector<int> overs{points[0].ID};
    //std::multiset<int> overs{points[0].ID};

    // for every adjacent transition point
    for (auto i = 1u; i < points.size(); ++i) 
    {
    	/*
    	std::cout << "Overs: ";
    	for (auto it = overs.begin() ; it != overs.end() ; ++it) std::cout << *it << " ";
    	std::cout << "\n";
    	*/


        auto &a = points[i - 1];
        auto &b = points[i];

        // if there is a jump in between transition points
		if (a.location < b.location)
		{
			switch (overs.size())
			{
				// no segment
				case 0 :
               		break;

				// ony one segment
               	case 1 : 
               		if (shouldNodeBeKeptIfItHasASegment)
               		{
               			if (this->newNodeID == -1)
							this->newNodeID = nodeTable.addNode(myNode);
               			B.push_back({a.location, b.location, this->newNodeID});	
               			edgeTable.addEdge({a.location, b.location, this->newNodeID, *overs.begin()});
               		} else
               		{
               			B.push_back({a.location, b.location, *overs.begin()});
               		}
               		break;
						

				// overlapping segment
				default :
					if (this->newNodeID == -1)
							this->newNodeID = nodeTable.addNode(myNode);
					
					for (auto it = overs.begin() ; it != overs.end() ; ++it)
					{
						assert(*it < this->newNodeID);
						edgeTable.addEdge({a.location, b.location, this->newNodeID, *it});
					}
					B.push_back({a.location, b.location, this->newNodeID});
					break;
			}
		}
        // update overlaps
        if (b.overlap)
			overs.insert(std::lower_bound(overs.begin(), overs.end(), b.ID), b.ID);
       		//overs.insert(b.ID);
        else
			overs.erase(std::lower_bound(overs.begin(), overs.end(), b.ID));
       		//overs.erase(overs.find(b.ID));
    }

	
    
    #ifdef DEBUG
    std::cout << "before merging adjacent:";
    segments = B;
    print();
    #endif
    // merge adjacent segments with same ID
	for (uint32_t i = 0 ; i < B.size()-1 ; ++i)
	{
		auto l = i;
		while (B.size() > l+1 && B[l].child == B[l+1].child && B[l].right == B[l+1].left)
		{
			++l;
		}
		if (l != i)
		{
			/*
			std::cout << "i = " << i << "\n";
			std::cout << "l= " << l<< "\n";
			*/
			auto newRight = B[l].right;
			assert(newRight > B[i].right);
			B[i].right = newRight;
			B.erase(B.begin() + i + 1, B.begin() + l + 1);
		}
	}

        
    segments.swap(B);
}

void T4TreeRec::HaplotypeOfSegments::reset_alreadySegmentized(bool b)
{
	alreadySegmentized = b;
}

void T4TreeRec::HaplotypeOfSegments::set_newNodeID(int id)
{
	newNodeID = id;
}


void T4TreeRec::HaplotypeOfSegments::simpleTransferSegmentsFromChild(T4TreeRec::HaplotypeOfSegments& childSegments, const T4TreeRec::Edge& edge, bool isOrdered)
{
	assert(childSegments.alreadySegmentized);


	//std::cout << "\tchildSegments.size() = " << childSegments.size() << "\n";
	for (uint32_t childSegmentIndex = 0 ; childSegmentIndex < childSegments.size() ; ++childSegmentIndex)
	{
		
		auto& childSegment = childSegments[childSegmentIndex];
		//std::cout << "\tCS: " << childSegment.left << " " << childSegment.right << "\n";

		if (edge.left  >= childSegment.right) continue;  // nothing from this childSegment will be used
		if (isOrdered && edge.right <= childSegment.left) break;

		auto movingSegment = getMovingSegment(edge, childSegment);

		//std::cout << "\tMS: {" << movingSegment.left << " " << movingSegment.right << "}\n";

		childSegments.removeMovingSequence(movingSegment, childSegmentIndex);

		segments.push_back(movingSegment); // Out of order is fine. It will be ordered after being segmentized	
		
	}
}


template<typename int_type> T4TreeRec::Segment& T4TreeRec::HaplotypeOfSegments::operator[](int_type& i)
{
	return segments[i];
}

uint32_t T4TreeRec::HaplotypeOfSegments::size()const{return segments.size();}
int T4TreeRec::HaplotypeOfSegments::getNewNodeID()const{return newNodeID;}
int T4TreeRec::HaplotypeOfSegments::getGeneration()const{return myNode.generation;}
int T4TreeRec::HaplotypeOfSegments::getPatchIndex() const{return myNode.patch;}
bool T4TreeRec::HaplotypeOfSegments::getAlreadySegmentized()const{return alreadySegmentized;}


T4TreeRec::HaplotypeOfSegments::HaplotypeOfSegments(std::vector<T4TreeRec::Segment>& segs, Node&& n)
:segments(segs), myNode(n), newNodeID(-1), alreadySegmentized(false)
{}
T4TreeRec::HaplotypeOfSegments::HaplotypeOfSegments(std::vector<T4TreeRec::Segment>& segs, Node& n)
:segments(segs), myNode(n), newNodeID(-1), alreadySegmentized(false)
{}
T4TreeRec::HaplotypeOfSegments::HaplotypeOfSegments(std::vector<T4TreeRec::Segment> segs, Node n)
:segments(segs), myNode(n), newNodeID(-1), alreadySegmentized(false)
{}
T4TreeRec::HaplotypeOfSegments::HaplotypeOfSegments(Node n)
:myNode(n), newNodeID(-1), alreadySegmentized(false)
{}

std::vector<T4TreeRec::Segment>& T4TreeRec::HaplotypeOfSegments::getSegmentsRef()
{
	return segments;
}

void T4TreeRec::HaplotypeOfSegments::transferSegmentsFromChild(T4TreeRec::HaplotypeOfSegments& childSegments, const T4TreeRec::Edge& edge, T4TreeRec::EdgeTable& newEdgeTable, T4TreeRec::NodeTable& newNodeTable, bool shouldNodeBeKeptIfItHasASegment)
{
	assert(edge.left < edge.right);
	assert(this->alreadySegmentized == false);

	/*
	#ifdef DEBUG
	std::cout << "Before segmentation\n";
	childSegments.print();
	#endif
	*/
	//std::cout << "edge: {" << edge.left << " " << edge.right << "}\n";



	if (!childSegments.alreadySegmentized) childSegments.segmentize(newNodeTable, newEdgeTable, shouldNodeBeKeptIfItHasASegment);	
	/*
	#ifdef DEBUG
	std::cout << "After segmentation\n";
	childSegments.print();
	#endif
	*/

	this->simpleTransferSegmentsFromChild(childSegments, edge);

	/*
	#ifdef DEBUG
	std::cout << "After transfer\n";
	childSegments.print();
	#endif
	*/
}


