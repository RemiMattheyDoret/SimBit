template<typename T>
class MemoryPool
{
	std::deque<T> objects; // because a deque never invalidates a pointer
	std::deque<T*> small_available;
	std::deque<T*> big_available;

	float multiplier = 1.1; // default multiplier value
	size_t nbNEWs = 0;
	size_t nbDELETESs = 0;
 
private:
	void reallocate();
	void increaseSize(size_t newSize);

public:
	MemoryPool();
	~MemoryPool();

	void resize(size_t newSize);
	size_t size() const;

	void makeBig(T* o);

	void setMultiplier(float m);
	void shrink_all_to_fit();

	T* NEW(bool isSmall);
	void DELETE(T* ptr, bool isSmall);

	void moveBigToSmall();
	void emptySmallMemory();

	void computeSomeStats();
};
