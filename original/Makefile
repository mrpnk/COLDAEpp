FORT = gfortran
FORTFLAGS = -Wfatal-errors -O3 -ffree-line-length-none -std=legacy

SRC := ./linpack/
OBJ := ./bin/

SOURCES := $(wildcard $(SRC)/*.f)
OBJECTS := $(patsubst $(SRC)/%.f, $(OBJ)%.o, $(SOURCES))

SOURCES90 := $(wildcard $(SRC)/*.f90)
OBJECTS90 := $(patsubst $(SRC)/%.f90, $(OBJ)%.o, $(SOURCES90))

all: test1 test2 test3 test4

test1: test1.f coldae.o $(OBJECTS) $(OBJECTS90)
	$(FORT) $(FORTFLAGS) $^ -o test1.exe
test2: test2.f coldae.o $(OBJECTS) $(OBJECTS90)
	$(FORT) $(FORTFLAGS) $^ -o test2.exe
test3: test3.f coldae.o $(OBJECTS) $(OBJECTS90)
	$(FORT) $(FORTFLAGS) $^ -o test3.exe
test4: test4.f90 colnew.o $(OBJECTS) $(OBJECTS90)
	$(FORT) $(FORTFLAGS) $^ -o test4.exe

bench1: bench1.f90 coldae.o $(OBJECTS) $(OBJECTS90)
	$(FORT) $(FORTFLAGS) $^ -o bench1.exe

coldae.o: coldae.f
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@
colsys.o: colsys.f
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@
colnew.o: colnew.f
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@

$(OBJ)%.o: $(SRC)/%.f
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@
	
$(OBJ)%.o: $(SRC)/%.f90
	$(FORT) $(FORTFLAGS) $< -I $(SRC) -c -o $@

clean:
	rm -f $(OBJ)*.o $(OBJECTS)


#     $@ means "name of the target"