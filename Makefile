SimBit: 
	g++ -std=c++14 src/main.cpp -o bin/SimBit -O3

noOpt: 
	g++ -std=c++14 src/main.cpp -o bin/SimBit -Wall

debug:
	g++ -std=c++14 src/main.cpp -o bin/SimBit -Wall -DDEBUG

debugSani:
	g++ -std=c++14 src/main.cpp -o bin/SimBit -Wall -DDEBUG -fsanitize=address

clean:
	rm bin/*
