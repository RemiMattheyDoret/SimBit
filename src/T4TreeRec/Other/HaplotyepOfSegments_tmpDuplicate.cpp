
class HaplotypeOfSegments
{
private:
	std::deque<Segment> segments;
	int id;


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

public:

	size_t nbSegments(){return segments.size();}
	void   set_id(int ID){id=ID;}
	int get_id(){return id;}

	HaplotypeOfSegments(std::deque<Segment>& segs, int& ID):segments(segs), id(ID){}
	HaplotypeOfSegments(std::deque<Segment> segs, int ID):segments(segs), id(ID){}
	HaplotypeOfSegments(int& ID):segments(segs), id(ID){}


	void transferSegmentsFromChild(HaplotypeOfSegments& childSegments, const Edge& edge, EdgeTable& edgeTable, NodeTable& nodeTable)
	{
		assert(childSegments.size());
		assert(edge.left < edge.right);
		assert(edge.child != this->id);

		int newNodeID= -1;
		size_t parentSegmentIndex = 0; 
		bool isSegmentsEmpty = segments.size() == 0;
		for (size_t childSegmentIndex = 0 ; childSegmentIndex < childSegments.size() ; ++childSegmentIndex)
		{
			auto& childSegment = childSegments[childSegmentIndex];



			///////////////////////
			// get movingSegment //
			///////////////////////

			if (edge.left  >= childSegment.right) continue; // nothing from this childSegment will be used
			if (edge.right <= childSegment.left) break;     // no need to look at other child segments

			auto& MSleft = std::max(edge.left, childSegment.left);
			auto& MSright = std::min(edge.right, childSegment.right);

			#ifdef DEBUG
			assert(MSleft < MSright);
			#endif

			movingSegment(
				MSleft,
				MSright,
				childSegment.child
			);





			/////////////////////////////////////
			// remove movingSegment from child //
			/////////////////////////////////////

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



			/*
			# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
			# //////////////////////////////////////////////////////// #
			# // Move the sequence and record in the tree if needed // #
			# //////////////////////////////////////////////////////// # 
			# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
			*/
	

			

			if (isSegmentsEmpty)
			{
				/////////////////////////////////////////////////////
				// Simple case where this->segments is still empty //
				/////////////////////////////////////////////////////

				segments.push_back(movingSegment);


			} else
			{
				////////////////////////////////////////////////////
				// Complex case where this->segments is not empty //
				////////////////////////////////////////////////////

				// Insert movingSegment up until there is nothing left in this movingSegment
				while (movingSegment.left < movingSegment.right)
				{
					/*
						Possibilities:
							currentParentSegment:    -----
							movingSegment:           -----

							currentParentSegment:    -----
							movingSegment:           --------  Watch out for the nextParentSegment

							currentParentSegment:    ----------
							movingSegment:           -------

							currentParentSegment:    -----
							movingSegment:             ---

							currentParentSegment:    -----
							movingSegment:             -----   Watch out for the nextParentSegment

							currentParentSegment:    ----------
							movingSegment:             -----


							currentParentSegment:    -----
							movingSegment:                  -----  Watch out for the nextParentSegment

							currentParentSegment:    -----
							movingSegment:                -----    Watch out for the nextParentSegment


						Note that neither currentParentSegment or nextParentSegment can start before movingSegment starts
					*/


					// Search previous parent sequence
					while(
						parentSegmentIndex < this->segments.size()
						&&
						this->segments[parentSegmentIndex].left < movingSegment.left
					) ++parentSegmentIndex;
					--parentSegmentIndex;
					#ifdef DEBUG
					assert(parentSegmentIndex > 0 && parentSegmentIndex < this->segments.size());
					assert(segments[parentSegmentIndex].left < movingSegment.left);
					#endif

					auto& currentParentSegment = this->segments[parentSegmentIndex];
					
					if (movingSegment.left >= currentParentSegment.right)
					{
						/*

							currentParentSegment:    -----
							nextParentSegment:               ??????????
							movingSegment:                  ------


							currentParentSegment:    -----
							nextParentSegment:            ??????????
							movingSegment:                ------

							currentParentSegment:    -----
							nextParentSegment:              ??
							movingSegment:                ------
						*/

						bool isIntersection = segments.size() > currentParentSegment + 1 && segments[currentParentSegment + 1].left < movingSegment.right;
						if (isIntersection)
						{
							if (newNodeID == -1) newNodeID = nodeTable.addNode(GP->CurrentGeneration);
							edgeTable.addEdge(segments[currentParentSegment + 1].left, segments[currentParentSegment + 1].right, segments[currentParentSegment + 1].child, newNodeID);
							edgeTable.addEdge(movingSegment, newNodeID);

							#ifdef DEBUG
							assert(segments[currentParentSegment + 1].left > movingSegment.left);
							#endif
							segments[currentParentSegment + 1].left = movingSegment.left;
						} else
						{
							++parentSegmentIndex;
							segments.insert(
								segments.begin() + nextParentSegment,
								{
									movingSegment.left,
									movingSegment.right,
									movingSegment.child
								}
							);
						}

					} else
					{
						/*
							currentParentSegment:    -----
							movingSegment:           -----

							currentParentSegment:    -----
							movingSegment:           --------

							currentParentSegment:    ----------
							movingSegment:           -------

							currentParentSegment:    -----
							movingSegment:             ---

							currentParentSegment:    -----
							movingSegment:             -----

							currentParentSegment:    ----------
							movingSegment:             -----
						*/


						if (movingSegment.left == currentParentSegment.left)
						{
							/*
								currentParentSegment:    -----
								movingSegment:           -----

								currentParentSegment:    -----
								movingSegment:           --------

								currentParentSegment:    ----------
								movingSegment:           -------
							*/

							if (movingSegment.right > currentParentSegment.right)
							{
								/*
									currentParentSegment:    -----
									movingSegment:           --------
								*/

								currentParentSegment.child = this->id;
								if (newNodeID == -1)
								{
									newNodeID = nodeTable.addNode(GP->CurrentGeneration);
								}
								edgeTable.addEdge({currentParentSegment.left, currentParentSegment.right, this->id, movingSegment.child});
								edgeTable.addEdge({currentParentSegment.left, currentParentSegment.right, this->id, currentParentSegment.child});

								movingSegment.left = currentParentSegment.right;
							} else if (movingSegment.right == currentParentSegment.right)
							{
								/*
									currentParentSegment:    -----
									movingSegment:           -----
								*/
								currentParentSegment.child = this->id;
								if (newNodeID == -1)
								{
									newNodeID = nodeTable.addNode(GP->CurrentGeneration);
									edgeTable.addEdge({movingSegment.left, movingSegment.right, this->id, movingSegment.child});
								}

							} else
							{
								/*
									currentParentSegment:    ----------
									movingSegment:           -------
								*/
								if (currentParentSegment.child != movingSegment.child)
								{
									// first insert the next element
									transferSegmentsFromChild_eventualInsertion(
										parentSegmentIndex, 
										{
											movingSegment.right,
											currentParentSegment.right,
											currentParentSegment.child
										},
										false
									);

									// Now change the current
									currentParentSegment.right = movingSegment.right;
									currentParentSegment.child = this->id;
									if (newNodeID == -1)
										newNodeID = nodeTable.addNode(GP->CurrentGeneration);
								}
								
								movingSegment.right = movingSegment.left; // set it to zero of length
							}

						} else
						{
							#ifdef DEBUG
							assert(movingSegment.left > currentParentSegment.left);
							#endif

							/*
								currentParentSegment:    -----
								movingSegment:             ---

								currentParentSegment:    -----
								movingSegment:              --------

								currentParentSegment:    ----------
								movingSegment:             -----
							*/

							if ( movingSegment.right > currentParentSegment.right )
							{
								/*
									currentParentSegment:    -----
									movingSegment:              --------
								*/

								if ( currentParentSegment.child != movingSegment.child )
								{
									transferSegmentsFromChild_eventualInsertion(
										parentSegmentIndex,
										{
											movingSegment.left,
											currentParentSegment.right,
											this->id
										},
										true
									);
									if (newNodeID == -1)
										newNodeID = nodeTable.addNode(GP->CurrentGeneration);

									currentParentSegment.right = movingSegment.left;
								}
								movingSegment.left = currentParentSegment.right;

							} else if ( movingSegment.right == currentParentSegment.right )
							{
								/*
									currentParentSegment:    -----
									movingSegment:             ---
								*/
								if ( currentParentSegment.child != movingSegment.child )
								{
									transferSegmentsFromChild_eventualInsertion(
										parentSegmentIndex,
										{
											movingSegment.left,
											movingSegment.right,
											this->id
										},
										true
									);
									if (newNodeID == -1)
										newNodeID = nodeTable.addNode(GP->CurrentGeneration);

									currentParentSegment.right = movingSegment.left;
								}
							} else
							{
								/*
									currentParentSegment:    ----------
									movingSegment:             -----
								*/
								if ( currentParentSegment.child != movingSegment.child )
								{
									transferSegmentsFromChild_eventualInsertion(
										parentSegmentIndex,
										{
											movingSegment.left,
											movingSegment.right,
											this->id
										},
										false
									);
									
									transferSegmentsFromChild_eventualInsertion(
										parentSegmentIndex,
										{
											movingSegment.right,
											currentParentSegment.right,
											this->id
										},
										true
									);
									if (newNodeID == -1)
										newNodeID = nodeTable.addNode(GP->CurrentGeneration);

									currentParentSegment.right = movingSegment.left;
									movingSegment.right = movingSegment.left;
								}
							}
						}

						
						#ifdef DEBUG
						assert(currentParentSegment.left < currentParentSegment.right);
						#endif
					}



					// Increment parentSegmentIndex
					++parentSegmentIndex // This should be optional though
				}
			} // End of if-else segments is empty
		} // end of for loop through childSegments

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