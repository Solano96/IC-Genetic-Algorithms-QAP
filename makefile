
CXX=g++
CXXFLAGS= -O2

main : main.cpp
	 $(CXX) $(CXXFLAGS) -o main main.cpp

clean:
	rm main
