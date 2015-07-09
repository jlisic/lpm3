# simple make file
#
CC=gcc
#CFLAGS=-g -DDDEBUG
CFLAGS=-g
LDFLAGS=-lm

lpm3: lpm3.o kdtree_lpm.o
	$(CC) -o lpm3 $(LDFLAGS) lpm3.o kdtree_lpm.o 

lpm3.o: lpm3.c
	$(CC) $(CFLAGS) -DCLI -c lpm3.c 

kdtree_lpm.o: kdtree_lpm.c
	$(CC) $(CFLAGS) -DCLI -c kdtree_lpm.c 

clean:
	rm -f lpm3 *.o *.so *.dll *.dylib
