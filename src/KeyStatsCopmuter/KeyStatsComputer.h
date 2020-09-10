class KeyStatsComputer
{
private:
	std::pair<std::vector<double>, std::vector<double>> compute_iHS(PopGenData& data);
	std::array<std::vector<double>, 3> compute_Hxy(PopGenData& data)
public:
	void computeAndWrite(PopGenData& data, std::string& outFilename);
};

