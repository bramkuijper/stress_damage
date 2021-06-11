


EXE=stress_damage.exe
CPP=stress_damage.cpp

EXE_LH=stress_damage_lh.exe
CPP_LH=stress_damage_lh.cpp

EXE_TEST=test.exe
CPP_TEST=test.cpp


CXX=g++
CXXFLAGS=-Wall -O3 -ggdb


all : $(EXE_TEST) $(EXE) $(EXE_LH) 

$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

$(EXE_LH) : $(CPP_LH)
	$(CXX) $(CXXFLAGS) -o $(EXE_LH) $(CPP_LH)


$(EXE_TEST) : $(CPP_TEST)
	$(CXX) $(CXXFLAGS) -o $(EXE_TEST) $(CPP_TEST)
