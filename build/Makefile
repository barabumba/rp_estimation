CC=icc -qopenmp
SOURCES=main.c processing.c psd.c sample_gen.c misc.c
HEADERS=defs.h processing.h psd.h sample_gen.h misc.h

CFLAGS=-std=c99 -c -O3 -Wall -I/home/turbin_mm/.gsl/include
LDFLAGS=-lm -lgsl -lgslcblas -L/home/turbin_mm/.gsl/lib

SRC_DIR = ../src/

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=rp_estimation

CSOURCES=$(addprefix $(SRC_DIR), $(SOURCES))

all: $(CSOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: $(SRC_DIR)%.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o

.PHONY: clean
