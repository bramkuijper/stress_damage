# makefile to compile the c++ files into executables
# Bram Kuijper 2021


EXE=stress_damage.exe
CPP=stress_damage.cpp

# life history model
EXE_LH=stress_damage_lh.exe
CPP_LH=stress_damage_lh.cpp

# terminal rewards
EXE_TERM=stress_damage_terminal.exe
CPP_TERM=stress_damage_terminal.cpp

# compiler is gcc
CXX=g++
CXXFLAGS=-Wall -O3 

all : $(EXE) $(EXE_LH) $(EXE_TERM)

$(EXE_TERM) : $(CPP_TERM)
	$(CXX) $(CXXFLAGS) -o $(EXE_TERM) $(CPP_TERM)

$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

$(EXE_LH) : $(CPP_LH)
	$(CXX) $(CXXFLAGS) -o $(EXE_LH) $(CPP_LH)

clean :
	rm -rf $(EXE)
	rm -rf $(EXE_LH)
	rm -rf $(EXE_TERM)
