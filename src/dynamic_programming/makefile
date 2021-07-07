


EXE=stress_damage.exe
CPP=stress_damage.cpp

EXE_LH=stress_damage_lh.exe
CPP_LH=stress_damage_lh.cpp

CXX=g++
CXXFLAGS=-Wall -O3 

all : $(EXE) $(EXE_LH) 

$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

$(EXE_LH) : $(CPP_LH)
	$(CXX) $(CXXFLAGS) -o $(EXE_LH) $(CPP_LH)


