FORT = gfortran
FORTFLAGS = -Wfatal-errors -w -std=legacy -O2 -Wall


SRC := ./
OBJ := ./

SOURCES := $(wildcard $(SRC)/*.f90)
OBJECTS := $(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SOURCES))

all: $(OBJECTS)
	$(FORT) $(FORTFLAGS) $^ -o f90test.exe
	
	
$(OBJ)/%.o: $(SRC)/%.f90
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@
	
