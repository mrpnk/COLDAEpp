FORT = gfortran
FORTFLAGS = -Wfatal-errors -std=legacy -O3 -Wall -std=legacy


SRC := ./
OBJ := ./

SOURCES := $(wildcard $(SRC)/*.f)
OBJECTS := $(patsubst $(SRC)/%.f, $(OBJ)/%.o, $(SOURCES))

SOURCES90 := $(wildcard $(SRC)/*.f90)
OBJECTS90 := $(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SOURCES90))

all: $(OBJECTS) $(OBJECTS90)
	$(FORT) $(FORTFLAGS) $^ -o ftest.exe
	
	
$(OBJ)/%.o: $(SRC)/%.f
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@
	
$(OBJ)/%.o: $(SRC)/%.f90
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@
	
clean:
	rm -f *.o $(OBJECTS)


#     $@ means "name of the target"