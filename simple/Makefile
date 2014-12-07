# Makefile

CC       = icc
CFLAGS   = -Wall -O3 -xHost -mkl
LDFLAGS  = $(CFLAGS)
LDLIBS   = 
INCLUDES = 
SRC      = hh.c
OBJS     = hh.o
PROG     = hh.out


all: $(PROG)

$(PROG): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDLIBS) $(LDFLAGS)

asm: $(INCLUDES) $(SRC)
	$(CC) -S $(SRC) $(CFLAGS)


*.o : $(INCLUDES) $(SRC)

.c.o:

clean:
	rm -fv $(OBJ) *~ *.o

