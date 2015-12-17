ARB_INCLUDES=-I. -I../arb -I/usr/local/include/flint
ARB_LIBS=-L../arb
ARB_LD_LIBRARY=LD_LIBRARY_PATH=/usr/local/bin:../arb

#VALGRIND=valgrind --leak-check=full --show-leak-kinds=all
VALGRIND=

CC=gcc
DEFS=
GLOBAL_DEFS=
WARNINGS=-Wall  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wmissing-declarations -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized

CFLAGS+= $(WARNINGS) -O3 $(DEFS) $(GLOBAL_DEFS) -march=native

# add a flag for valgrind
CFLAGS+= -g

TEST_EXECUTABLES=bin/t-factor_refinement \
		 bin/t-femtocas \
		 bin/t-expressions \
		 bin/t-generators

CHECKS=check_factor_refinement \
       check_femtocas \
       check_expressions \
       check_generators

all: factor_refinement.o femtocas.o expressions.o generators.o

# force the test scripts to be built by listing them as requirements
tests: $(TEST_EXECUTABLES)

# run the test scripts
check: $(CHECKS)


# factor refinement

check_factor_refinement: bin/t-factor_refinement
	valgrind bin/t-factor_refinement

factor_refinement.o:
	$(CC) factor_refinement.c -c $(CFLAGS) -lflint -lgmp

bin/t-factor_refinement: tests/t-factor_refinement.c factor_refinement.o
	$(CC) tests/t-factor_refinement.c factor_refinement.o \
		-o bin/t-factor_refinement \
		-I. $(CFLAGS) -lflint -lgmp


# femtocas

check_femtocas: bin/t-femtocas
	$(ARB_LD_LIBRARY) valgrind bin/t-femtocas

femtocas.o: femtocas.c
	$(CC) femtocas.c -c $(ARB_INCLUDES) $(CFLAGS) -lflint -lgmp -larb

bin/t-femtocas: tests/t-femtocas.c femtocas.o
	$(CC) tests/t-femtocas.c femtocas.o \
		-o bin/t-femtocas \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# expressions

check_expressions: bin/t-expressions
	$(ARB_LD_LIBRARY) valgrind bin/t-expressions

expressions.o: expressions.c
	$(CC) expressions.c -c $(ARB_INCLUDES) $(CFLAGS) femtocas.o \
		-lflint -lgmp -larb

bin/t-expressions: tests/t-expressions.c expressions.o
	$(CC) tests/t-expressions.c expressions.o femtocas.o \
		-o bin/t-expressions \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# generators

check_generators: bin/t-generators
	$(ARB_LD_LIBRARY) valgrind bin/t-generators

generators.o: generators.c
	$(CC) generators.c -c $(ARB_INCLUDES) $(CFLAGS) \
		femtocas.o expressions.o \
		-lflint -lgmp -larb

bin/t-generators: tests/t-generators.c generators.o
	$(CC) tests/t-generators.c generators.o femtocas.o expressions.o \
		-o bin/t-generators \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# double precision wavefront implementation

wavefront_double.o: wavefront_double.c
	$(CC) wavefront_double.c -c $(CFLAGS) \
		-lflint -lgmp


# integer vector wavefront implementation

wavefront_hermite.o: wavefront_hermite.c
	$(CC) wavefront_hermite.c -c $(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb


# this is the command line binary executable

bin/arbtkf91: arbtkf91.c \
	femtocas.o expressions.o generators.o \
	wavefront_double.o wavefront_hermite.o
	$(CC) arbtkf91.c \
		femtocas.o expressions.o generators.o \
		wavefront_double.o wavefront_hermite.o \
		-o bin/arbtkf91 \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb -lm


# run an example

example1a: bin/arbtkf91
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/arbtkf91 \
		--sequence-1 ACGATA \
		--sequence-2 AGTGGTA \
		--lambda-num 1 --lambda-den 1 \
		--mu-num 2 --mu-den 1 \
		--tau-num 1 --tau-den 10 \
		--pa-num 1 --pa-den 4 \
		--pc-num 1 --pc-den 4 \
		--pg-num 1 --pg-den 4 \
		--pt-num 1 --pt-den 4

example1b: bin/arbtkf91
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/arbtkf91 \
		--sequence-1 ACGATA \
		--sequence-2 AGTGGTA \
		--lambda-num 1 --lambda-den 1 \
		--mu-num 2 --mu-den 1 \
		--tau-num 1 --tau-den 10 \
		--pa-num 27 --pa-den 100 \
		--pc-num 24 --pc-den 100 \
		--pg-num 26 --pg-den 100 \
		--pt-num 23 --pt-den 100

example2a: bin/arbtkf91
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/arbtkf91 \
		--sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA \
		--sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA \
		--lambda-num 1 --lambda-den 1 \
		--mu-num 2 --mu-den 1 \
		--tau-num 1 --tau-den 10 \
		--pa-num 1 --pa-den 4 \
		--pc-num 1 --pc-den 4 \
		--pg-num 1 --pg-den 4 \
		--pt-num 1 --pt-den 4

example2b: bin/arbtkf91
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/arbtkf91 \
		--sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA \
		--sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA \
		--lambda-num 1 --lambda-den 1 \
		--mu-num 2 --mu-den 1 \
		--tau-num 1 --tau-den 10 \
		--pa-num 27 --pa-den 100 \
		--pc-num 24 --pc-den 100 \
		--pg-num 26 --pg-den 100 \
		--pt-num 23 --pt-den 100



clean:
	rm -f *.o
	rm -f bin/example
	rm -f bin/t-factor_refinement
	rm -f bin/t-femtocas
	rm -f bin/t-expressions
	rm -f bin/t-generators
