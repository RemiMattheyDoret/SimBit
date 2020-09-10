
class HaplotypeOfSegments
{
private:
	std::deque<Segment> segments;
	int oldNodeID;
	int newNodeID = -1;



	void assertOrderingAndNoOverlap()
	{
		for (size_t i = 0 ; i < segments.size() ; ++i)
		{
			auto& segment = segments[i];
			assert(segment.left < segment.right);
			if (i)
			{
				assert(segment.left >= segments[i].right);
			}
		}
	}

	static void removeMovingSequenceFromChildSegment(Segment& movingSegment, std::deque<Segment>& childSegments, size_t& childSegmentIndex)
	{
		if (movingSegment.left > childSegment.left && movingSegment.right < childSegment.right)
		{
			/*
				childSegment    ------------
				MSeq               -----

				Need to duplicate
			*/

			// First insert new childSegment
			childSegments.insert(
				childSegments.begin() + childSegmentIndex,
				{
					movingSegment.right,
					childSegments.right,
					childSegments.child
				}
			);
			#ifdef DEBUG
			assert(movingSegment.right < childSegments.right);
			#endif

			// Modify current childSegment
			childSegments.right = movingSegment.left;
			#ifdef DEBUG
			assert(childSegments.left < childSegments.right);
			#endif

			// increment childSegmentIndex so that we don't look at it again
			++childSegmentIndex;

			// That should work find because there should not be any overlap in childSegments
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

			if (childSegment.left < movingSegment.left)
			{
				/*
					childSegment    -------
					MSeq               ----
				*/
				#ifdef DEBUG
				assert(movingSegment.right == childSegments.right);
				#endif
				childSegments.right = movingSegment.left;
			} else if (childSegment.right > movingSegment.right)
			{
				/*
					childSegment       -------
					MSeq               ----
				*/
				#ifdef DEBUG
				assert(movingSegment.left == childSegment.left);
				#endif
				childSegment.left = movingSegment.right;
			} else
			{
				/*
					childSegment       -------
					MSeq               -------
				*/
				#ifdef DEBUG
				assert(movingSegment.left == childSegment.left);
				assert(movingSegment.right == childSegment.right);
				#endif

				childSegments.erase(childSegments.begin() + childSegmentIndex);
				--childSegmentIndex;
				// This is a a little dangerous because childSegment pointer can be non-sense now. But I just have to not use it from now on.
			}
		}
	}


	Segment getMovingSegment(Edge& edge, Segment& childSegment)
	{
		if (edge.left  >= childSegment.right) continue; // nothing from this childSegment will be used
			if (edge.right <= childSegment.left) break;     // no need to look at other child segments

			auto& MSleft = std::max(edge.left, childSegment.left);
			auto& MSright = std::min(edge.right, childSegment.right);

			#ifdef DEBUG
			assert(MSleft < MSright);
			#endif

			retturn Segment(
				MSleft,
				MSright,
				childSegment.child
			);
	}

public:

	size_t nbSegments(){return segments.size();}
	int getNewNodeID(){return newNodeID;}
	int getOldNodeID(){return oldNodeID;}

	HaplotypeOfSegments(std::deque<Segment>& segs, int& ID):segments(segs), id(ID){}
	HaplotypeOfSegments(std::deque<Segment> segs, int ID):segments(segs), id(ID){}
	HaplotypeOfSegments(int& ID):segments(segs), id(ID){}


	void transferSegmentsFromChild(HaplotypeOfSegments& childSegments, const Edge& edge, EdgeTable& edgeTable, NodeTable& nodeTable)
	{
		assert(childSegments.size());
		assert(edge.left < edge.right);
		assert(edge.child != this->oldNodeID);

		if (segments.size() == 0)
		{
			for (size_t childSegmentIndex = 0 ; childSegmentIndex < childSegments.size() ; ++childSegmentIndex)
			{
				auto& childSegment = childSegments[childSegmentIndex];

				auto movingSegment = getMovingSegment(edge, childSegment );

				HaplotypeOfSegments::removeMovingSequenceFromChildSegments(movingSegment, childSegments, childSegmentIndex);

				segments.push_back(movingSegment);
			}

			#ifdef DEBUG
			assert(!isTheSecondEdgeReceived);
			#endif
			isTheSecondEdgeReceived = true;
		} else
		{
			std::vector<Segment> C2;
			for (auto& childSegment : childSegments) C2.push_back(getMovingSegment(edge, childSegment ););

			std::vector<Segment> C1;
			C1.swap(segments);
			size_t C1i = 0;
			size_t C2i = 0;
			#ifdef DEBUG
			assert(C1.size());
			assert(C2.size());
			#endif
			
			segments.reserve(C1.size() + C2.size());

			Segment defaultMergedSegment(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::min(), -1);
			Segment mergedSegment = defaultMergedSegment;


			bool currentC1IntersectsWithNextC2 = false;
			bool currentC2IntersectsWithNextC1 = false;

			while ( true )
			{
		        if (C2i == C2.size() || C1i == C1.size())
		        {
		            break;
		        }

		        if (C2[C2i] < C1[C1i])
		        {
		        	if (C2[C2i].right < C1[C1i].left)
		        	{
		        		// intersection with the next one
		        		// set left of mergedSegment
		        		if (newNodeID == -1) newNodeID = nodeTable.addNode(GP->CurrentGeneration);
		        		if (mergedSegment.child == -1)
		        		{
		        			mergedSegment.child = this->newNodeID;
		        			mergedSegment.left  = C2[C2i].left;
		        		} else
		        		{
		        			#ifdef DEGUG
		        			assert(mergedSegment.child == this->newNodeID);
		        			assert(mergedSegment.left <= C2[C2i].left);
		        			#endif
		        		}
		        		currentC2intersect = true;
		        		
		        	} else
		        	{
		        		if (mergedSegment.child == -1)
		        		{
		        			segments.push_back(C2[C2i]);
		        		} else
		        		{
		        			// intersection with the previous one but not with the next one
		        			#ifdef DEBUG
		        			assert(this->newNodeID != -1);
		        			assert(mergedSegment.child == this->newNodeID);
		        			#endif
		        			mergedSegment.right = C2[C2i].right;
		        			segments.push_back(mergedSegment);
		        			mergedSegment = defaultMergedSegment;
		        		}
		        	}
		        	++C2i;
		        } else
		        {
		        	if (C1[C1i].right < C2[C2i].left)
		        	{
		        		// intersection with the next one
		        		// set left of mergedSegment
		        		if (newNodeID == -1) newNodeID = nodeTable.addNode(GP->CurrentGeneration);
		        		if (mergedSegment.child == -1)
		        		{
		        			mergedSegment.child = this->newNodeID;
		        			mergedSegment.left  = C1[C1i].left;
		        		} else
		        		{
		        			#ifdef DEGUG
		        			assert(mergedSegment.child == this->newNodeID);
		        			assert(mergedSegment.left <= C1[C1i].left);
		        			#endif
		        		}
		        		currentC1intersect = true;
		        		
		        	} else
		        	{
		        		if (mergedSegment.child == -1)
		        		{
		        			segments.push_back(C1[C1i]);
		        		} else
		        		{
		        			// intersection with the previous one but not with the next one
		        			#ifdef DEBUG
		        			assert(this->newNodeID != -1);
		        			assert(mergedSegment.child == this->newNodeID);
		        			#endif
		        			mergedSegment.right = C1[C1i].right;
		        			segments.push_back(mergedSegment);
		        			mergedSegment = defaultMergedSegment;
		        		}
		        	}
		        	++C1i;
		        }

		        // Build table
		        #ifdef DEBUG
		        assert(!(currentC1intersect && currentC2intersect))
		        #endif		        	

		        if (currentC1intersect)
		        {
		        	#ifdef DEBUG  
		        	assert(newNodeID != -1); 
		        	#endif
		        	if (isTheSecondEdgeReceived) edgeTable.addEdge(C1[C1i], newNodeID);
		        	currentC1intersect = false;
		        } else if (currentC2intersect)
		        {
		        	#ifdef DEBUG  
		        	assert(newNodeID != -1); 
		        	#endif
		        	edgeTable.addEdge(C2[C2i], newNodeID);
		        	currentC2intersect = false;
		        }
		    }

		    #ifdef
		    assert(C2i == C2.size() || C1i == C1.size());
		    assert(C2i != C2.size() || C1i != C1.size());
		    #endif


		    // Copy what's left and complete mergedSegment if it was found to intersct with the current

		    if (mergedSegment.child != -1 && newNodeID == -1)
		    {
		    	newNodeID = nodeTable.addNode(GP->CurrentGeneration);
		    }

	    	if (C1i == C1.size())
		    {
		    	if (mergedSegment.child != -1)
		    	{
		    		mergedSegment.right = C2[C2i].right;
			    	segments.push_back(mergedSegment);

			    	if (isTheSecondEdgeReceived) EdgeTable.addEdge(C2[C2i], newNodeID);
		    	}
		    	std::copy(C2.begin() + C2i, C2.end(), segments.end());
		    } else
		    {
		    	if (mergedSegment.child != -1)
		    	{
		    		mergedSegment.right = C1[C1i].right;
			    	segments.push_back(mergedSegment);

			    	EdgeTable.addEdge(C1[C1i], newNodeID);
		    	}
		    	std::copy(C1.begin() + C1i, C1.end(), segments.end());
		    }


			isTheSecondEdgeReceived = false;
	    }

		#ifdef DEBUG
		assertOrderingAndNoOverlap();
		#endif
	}
};



/*
class SegmentContainer
{
	std::priority_queue<Segment, std::vector<Segment>, std::greater<Segment>> data;

	const Segment& top()
	{
		return data.top();
	}

	void pop()
	{
		data.pop();
	}

	void push(Segment& seg)
	{
		data.push(seg);
	}

	void push(Segment& seg)
	{
		data.push(seg);
	}
};
*/