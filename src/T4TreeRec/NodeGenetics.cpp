
/*T4TreeRec::NodeGenetics::NodeGenetics(){}
T4TreeRec::NodeGenetics::NodeGenetics(int g):birth(g){}
T4TreeRec::NodeGenetics::NodeGenetics(const NodeGenetics& n):mutations(n.mutations), birth(n.birth){}
T4TreeRec::NodeGenetics::NodeGenetics(const NodeGenetics&& n):mutations(n.mutations), birth(n.birth){}
T4TreeRec::NodeGenetics::NodeGenetics(std::vector<uint32_t> m, int b):mutations(m), birth(b){}
//void assertWasFullySet() const;
T4TreeRec::NodeGenetics T4TreeRec::NodeGenetics::operator=(const T4TreeRec::NodeGenetics& other)
{
	mutations = other.mutations;
	birth = other.birth;
	return *this;
}
*/

T4TreeRec::NodeGenetics::NodeGenetics(int id)
:birth(-1), ID(id)
{}

T4TreeRec::NodeGenetics::NodeGenetics()
:birth(-1), ID(-1)
{}

T4TreeRec::NodeGenetics::NodeGenetics(std::vector<uint32_t> m, int b, int id)
:mutations(m), birth(b), ID(id)
{}

T4TreeRec::NodeGenetics::NodeGenetics(std::vector<uint32_t>& m, int b, int id)
:mutations(m), birth(b), ID(id)
{}

int T4TreeRec::NodeGenetics::getID() const
{
	return ID;
}

void T4TreeRec::NodeGenetics::resetID(int newID)
{
	ID = newID;
}

int T4TreeRec::NodeGenetics::getGeneration() const
{
	return birth;
}

void T4TreeRec::NodeGenetics::setGeneration(const int g)
{
    birth = g;
}

void T4TreeRec::NodeGenetics::clear()
{
	mutations.clear();
	//whatHasBeenSetYet = std::priority_queue<Segment, std::vector<Segment>, std::less<Segment>>();
}

void T4TreeRec::NodeGenetics::shrink_to_fit()
{
	mutations.shrink_to_fit();
}

std::vector<uint32_t>& T4TreeRec::NodeGenetics::getMutations()
{
	return mutations;
}

const std::vector<uint32_t>& T4TreeRec::NodeGenetics::getMutations() const
{
    return mutations;
}

std::vector<uint32_t>&& T4TreeRec::NodeGenetics::moveMutations()
{
	return std::move(mutations);
}

template<typename INT> T4TreeRec::NodeGenetics T4TreeRec::NodeGenetics::getSubsetMutations(INT from, INT to) const
{
    assert(from < to);
    auto l = std::lower_bound(mutations.begin(), mutations.end(), from);
    auto r = std::lower_bound(l, mutations.end(), to);
    return NodeGenetics({l,r}, this->getGeneration(), this->getID());
}

void T4TreeRec::NodeGenetics::push_backMutations(std::vector<uint32_t>& m)
{
    if (mutations.size() && m.size()) assert(mutations.back() < m.front());
    mutations.insert(mutations.end(), m.begin(), m.end());
}

void T4TreeRec::NodeGenetics::push_backMutations(std::vector<uint32_t>&& m)
{
    if (mutations.size() && m.size()) assert(mutations.back() < m.front());
    mutations.insert(mutations.end(), m.begin(), m.end());
}

template<typename INT>
void T4TreeRec::NodeGenetics::mutateLocus(INT MutPosition)
{
	auto position = std::lower_bound(mutations.begin(), mutations.end(), MutPosition);

    if (position == mutations.end())
    {
        // not found
        if (SSP->T4_mutDirection != 0)
            mutations.push_back(MutPosition);  
    } else
    {
        if ( MutPosition == (*position))
        {
            // found
            if (SSP->T4_mutDirection != 1)
                mutations.erase(position);
        } else
        {
            // not found
            if (SSP->T4_mutDirection != 0)
                mutations.insert(position, MutPosition);
        }
    }
}

