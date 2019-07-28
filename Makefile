CC=gcc -fopenmp  -I/home/turbin_mm/.gsl/include -L/home/turbin_mm/.gsl/lib
CFLAGS=-std=c99 -c -Wall 
LDFLAGS=-lm -lgsl -lgslcblas 
SOURCES=main.c processing.c psd.c sample_gen.c
HEADERS=defs.h processing.h psd.h sample_gen.h
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=rp_estimation

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o: 
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm *.o

.PHONY: clean
