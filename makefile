CC=g++
CFLAGS=-I.
HEADERS = Header.h arrays.h functions.h optimize.h
OBJ =  main.cpp functions.cpp optimize.cpp fourier.cpp arrays.cpp plotandprint.cpp

CalcFarField: $(OBJ) $(HEADERS)
	$(CC) -o CalcFarField $(OBJ) $(CFLAGS)

