# Makefile

CC       = gcc
CFLAGS   = -Wall -O2
LDFLAGS  = $(CFLAGS)
LDLIBS   = -lm
INCLUDES = 
SRC      = hh.c
OBJS     = hh.o
PROG     = hh.out


all: $(PROG)

$(PROG): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDLIBS) $(LDFLAGS)

*.o : $(INCLUDES) $(SRC)

.c.o:

clean:
	rm -fv $(OBJ) *~ *.o

