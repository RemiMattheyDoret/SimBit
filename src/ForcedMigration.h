class ForcedMigration
{
private:
	std::vector<size_t> earlytimes;
	std::vector<bool> isCopy_atEarlyTimes;
	std::vector<size_t> from_atEarlyTimes;
	std::vector<size_t> to_atEarlyTimes;
	std::vector<size_t> nbInds_atEarlyTimes;

	std::vector<size_t> latetimes;
	std::vector<bool> isCopy_atLateTimes;
	std::vector<size_t> from_atLateTimes;
	std::vector<size_t> to_atLateTimes;
	std::vector<size_t> nbInds_atLateTimes;

public:
	void readInput(InputReader& input);	
	void forceMigrationIfNeeded(Pop& pop, bool isEarly);
};