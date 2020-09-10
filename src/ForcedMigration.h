class ForcedMigration
{
	std::vector<size_t> times;
	std::vector<std::vector<std::vector<size_t>>> froms_atTimes;
	std::vector<std::vector<std::vector<size_t>>> tos_atTimes;
	std::vector<std::vector<std::vector<size_t>>> nbHaploss_atTimes;

	void readInput(InputReader& input);
}