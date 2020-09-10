void T4TreeRec::Edge::print() const
{
	std::cout << "{" << left << " " << right << " " <<  parent << " "  << child << "}\n";
}

void T4TreeRec::Edge::swap(T4TreeRec::Edge& other)
{
	std::swap(left, other.left);
	std::swap(right, other.right);
	std::swap(parent, other.parent);
	std::swap(child, other.child);
}


bool T4TreeRec::Edge::operator>(const T4TreeRec::Edge& other) const
{
	if (this->parent == other.parent)
	{
		if (this->child == other.child)
		{
			return this->left < other.left;
		} else
		{
			return this->child > other.child;
		}
	} else
	{
		return this->parent > other.parent;
	}
}

bool T4TreeRec::Edge::operator<(const T4TreeRec::Edge& other) const
{
	if (this->parent == other.parent)
	{
		if (this->child == other.child)
		{
			return this->left > other.left;
		} else
		{
			return this->child < other.child;
		}
	} else
	{
		return this->parent < other.parent;
	}
}
