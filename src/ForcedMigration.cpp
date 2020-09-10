
void ForcedMigration::readInput(InputReader& input)
{
	while (input.IsThereMoreToRead())
	{
		if (input.peekNextElementString() != "times")
		{
			std::cout << "For option --forcedMigration, expected keyword 'times' but received " << input.peekNextElementString() << " instead.\n";
			abort();
		}

		while (input.peekNextElementString() != "from")
		{
			auto t = input.getNextElementInt();
			if (t <= 0 || t > GP->NbGenerations)
			{
				std::cout << "For option --forcedMigration, received a time value that is either not strictly positive or greater than the total number of generations. Received " << t << ".\n";
				abort();
			}
			times.push_back(t);
		}

		auto from = input.getNextElementInt();
		auto to = input.getNextElementInt();
	}
}

