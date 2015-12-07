LIBS=-lflint -lgmp


CC=gcc
DEFS=
GLOBAL_DEFS=
WARNINGS=-Wall  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wmissing-declarations -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized

CFLAGS+= $(WARNINGS) -O3 $(DEFS) $(GLOBAL_DEFS) -march=native

TARGET=factor_refinement.o


all: $(TARGET)

# force the test scripts to be built by listing them as requirements
tests: $(TARGET) tests/t-factor_refinement

# run the test scripts
check: tests/t-factor_refinement
	bin/t-factor_refinement

# build the test scripts
tests/t-factor_refinement: tests/t-factor_refinement.c $(TARGET)
	$(CC) tests/t-factor_refinement.c $(TARGET) \
		-o bin/t-factor_refinement \
		-I. $(CFLAGS) $(LIBS)

# build the factor refinement object file
$(TARGET):
	$(CC) factor_refinement.c -c $(CFLAGS) $(INCLUDES) $(LIBS)

clean:
	rm -f *.o
	rm -f $(TARGET)
	rm -f bin/t-factor_refinement
