CC = gcc

CCFLAGS = -Wall -O3

all: nbody-satellites analyze analyze-energy analyze-angular-momentum analyze-spinorbit-invariable

HEADERS = quaternion.h nbody.h 
CODE_LIBRARY = quaternion.c vector-stuff.c 
CODE = nbody-map.c nbody-io.c rigid.c nbody-energy.c nbody-angular-momentum.c jacobi.c corrector.c universal.c kepcart-new.c 

nbody-satellites: nbody-satellites.c $(HEADERS) $(CODE_LIBRARY) $(CODE)
	$(CC) $(CCFLAGS) -o nbody-satellites nbody-satellites.c -lm

analyze: analyze.c $(HEADERS) $(CODE_LIBRARY) $(CODE)
	$(CC) $(CCFLAGS) -o analyze analyze.c -lm

analyze-energy: analyze-energy.c $(HEADERS) $(CODE_LIBRARY) $(CODE)
	$(CC) $(CCFLAGS) -o analyze-energy analyze-energy.c -lm

analyze-angular-momentum: analyze-angular-momentum.c $(HEADERS) $(CODE_LIBRARY) $(CODE)
	$(CC) $(CCFLAGS) -o analyze-angular-momentum analyze-angular-momentum.c -lm

analyze-spinorbit-invariable: analyze-spinorbit-invariable.c $(HEADERS) $(CODE_LIBRARY) $(CODE)
	$(CC) $(CCFLAGS) -o analyze-spinorbit-invariable analyze-spinorbit-invariable.c -lm

clean:
	rm nbody-satellites analyze analyze-energy analyze-angular-momentum analyze-spinorbit-invariable



