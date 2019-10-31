class T5Mutation
{
	static unsigned long long currentID;

	unsigned long long ID;
	int generation;
	int locus; // That's a T5Locus

	T5Mutation(int locus);
	operator==(T5Mutation other);
	operator==(int otherLocus);
};
