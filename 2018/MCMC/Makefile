IDIR =../include
CC=gcc
CFLAGS=-v -lm -DL=16 -DHEATBATH -o Ising-L16

ODIR=obj
LDIR =../lib

DEPS = $(wildcard *.c)

OBJ = $(patsubst %.c,%.o,$(DEPS))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
