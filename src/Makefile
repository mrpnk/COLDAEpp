CPP = g++
CPPFLAGS = -std=gnu++2a -O2 -Wall
# -frounding-math -fsignaling-nans -ffast-math
#-Wfatal-errors

SRC := ./
OBJ := ./out/

SRC_LP := ./linpack/
OBJ_LP := ./out/

SOURCES := $(wildcard $(SRC)*.cpp) 
OBJECTS := $(patsubst $(SRC)%.cpp, $(OBJ)/%.o, $(SOURCES))

SOURCES_LP := $(wildcard $(SRC_LP)/*.cpp) 
OBJECTS_LP := $(patsubst $(SRC_LP)%.cpp, $(OBJ_LP)/%.o, $(SOURCES_LP))

all: $(OBJECTS) $(OBJECTS_LP) 
	$(CPP) $(CPPFLAGS) $^ -o cpptest.exe
	
	
$(OBJ)/%.o: $(SRC)/%.cpp
	$(CPP) $(CPPFLAGS) $< -I $(SRC) -c -o $@

$(OBJ_LP)/%.o: $(SRC_LP)/%.cpp
	$(CPP) $(CPPFLAGS) $< -I $(SRC_LP) -c -o $@
	
clean:
	rm -f *.o $(OBJECTS)


#     $@ means "name of the target"