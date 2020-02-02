
CXX=g++
CXXFLAGS= -O2

main : src/main.cpp
	 $(CXX) $(CXXFLAGS) -o main src/main.cpp

clean:
	rm main
