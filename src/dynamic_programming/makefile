# makefile to compile the c++ files into executables
# Bram Kuijper 2021


EXE=stress_damage.exe
CPP=stress_damage.cpp

EXE_LH=stress_damage_lh.exe
CPP_LH=stress_damage_lh.cpp

# compiler is gcc
CXX=g++
CXXFLAGS=-Wall -O3 

all : $(EXE) $(EXE_LH) 

$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

$(EXE_LH) : $(CPP_LH)
	$(CXX) $(CXXFLAGS) -o $(EXE_LH) $(CPP_LH)

clean :
	rm -rf $(EXE)
	rm -rf $(EXE_LH)
