### Remi Matthey-Doret
# You will need boost to compile SimBit. If the compilation fails with an error of the kind "src/main.cpp:55:38: fatal error: boost/algorithm/string.hpp: No such file or directory", please download the latest boost libraries and just place it in the SimBit folder
# You will need the standards of C++11 (or more recent) to compile SimBit.

# For the cluster orcinus.westgrid, you will have to do 
#      module load gcc/4.4.7/gsl-1.16
#      module load gcc/4.8.2rev203690
#      module load intel/14.0/boost-1.60.0
#      make c++11

SimBit: 
	g++ -std=c++14 src/main.cpp -o bin/SimBit -O3

c++11:
	g++ -std=c++11 src/main.cpp -o bin/SimBit -O3	

noOpt: 
	g++ -std=c++14 src/main.cpp -o bin/SimBit

debug:
	g++ -std=c++14 src/main.cpp -o bin/SimBit -DDEBUG

debugSani:
	g++ -std=c++14 src/main.cpp -o bin/SimBit -DDEBUG -fsanitize=address

sani:
	g++ -std=c++14 src/main.cpp -o bin/SimBit -fsanitize=address

clean:
	rm bin/*

