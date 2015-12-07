LIBS=-lflint -lgmp


CC=gcc
DEFS=
GLOBAL_DEFS=
WARNINGS=-Wall  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wmissing-declarations -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized

CFLAGS+= $(WARNINGS) -O3 $(DEFS) $(GLOBAL_DEFS) -march=native

TARGET=bin/my-executable

valgrind: CFLAGS+= -g



all:	$(TARGET)

valgrind: $(TARGET)
	valgrind --leak-check=yes $(TARGET)

$(TARGET):
	$(CC) factor_refinement.c -o $(TARGET) $(CFLAGS) $(INCLUDES) $(LIBS)

clean:
	rm -f *.o
	rm -f $(TARGET)
