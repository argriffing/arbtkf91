ARB_INCLUDES=-I. -I../arb -I/usr/local/include/flint
ARB_LD_LIBRARY=LD_LIBRARY_PATH=/usr/local/bin:../arb/build

CC=gcc
DEFS=
GLOBAL_DEFS=
WARNINGS=-Wall  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wmissing-declarations -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized

CFLAGS+= $(WARNINGS) -O3 $(DEFS) $(GLOBAL_DEFS) -march=native


all: factor_refinement.o femtocas.o

# force the test scripts to be built by listing them as requirements
tests: tests/t-factor_refinement tests/t-femtocas

# run the test scripts
check: check_factor_refinement check_femtocas


# factor refinement

check_factor_refinement: bin/t-factor_refinement
	bin/t-factor_refinement

factor_refinement.o:
	$(CC) factor_refinement.c -c $(CFLAGS) -lflint -lgmp

tests/t-factor_refinement: tests/t-factor_refinement.c factor_refinement.o
	$(CC) tests/t-factor_refinement.c factor_refinement.o \
		-o bin/t-factor_refinement \
		-I. $(CFLAGS) -lflint -lgmp


# femtocas

check_femtocas: bin/t-femtocas
	$(ARB_LD_LIBRARY) bin/t-femtocas

femtocas.o:
	$(CC) femtocas.c -c $(ARB_INCLUDES) $(CFLAGS) -lflint -lgmp

tests/t-femtocas: tests/t-femtocas.c femtocas.o
	$(CC) tests/t-femtocas.c femtocas.o \
		-o bin/t-femtocas \
		$(ARB_INCLUDES) $(CFLAGS) -lflint -lgmp


clean:
	rm -f *.o
	rm -f factor_refinement.o
	rm -f femtocas.o
	rm -f bin/t-factor_refinement
	rm -f bin/t-femtocas
