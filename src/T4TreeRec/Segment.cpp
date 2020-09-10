T4TreeRec::Segment::Segment(uint32_t l, uint32_t r, int c)
:left(l), right(r), child(c)
{}

void T4TreeRec::Segment::print() const
{
	std::cout << "{" << left << " " << right << " " << child << "}\n";
}