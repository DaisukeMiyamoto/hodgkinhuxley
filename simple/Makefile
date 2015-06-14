# Makefile

CC       = mpifccpx
CFLAGS   = -Kfast,openmp,optmsg=2 -DKCOMPUTER -Nrt_tune,src,sta
#CFLAGS   = -Kfast,openmp,optmsg=2 -DKCOMPUTER -Nrt_tune,src,sta -DUSE_FLOAT
LDFLAGS  = $(CFLAGS)
LDLIBS   = -lm
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

