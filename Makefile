ARB_INCLUDES=-I. -I../arb -I/usr/local/include/flint
ARB_LIBS=-L../arb
ARB_LD_LIBRARY=LD_LIBRARY_PATH=/usr/local/bin:../arb

VALGRIND=valgrind --leak-check=full --show-leak-kinds=all
#VALGRIND=

CC=gcc
DEFS=
GLOBAL_DEFS=
WARNINGS=-Wall  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wmissing-declarations -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized

#CC=g++
#DEFS=
#GLOBAL_DEFS=
#WARNINGS=-W -Wall -Wextra -pedantic -fpermissive


CFLAGS+= $(WARNINGS) -O3 $(DEFS) $(GLOBAL_DEFS) -march=native
CFLAGS+= -g
#CFLAGS+= -std=c++11


TEST_EXECUTABLES=bin/t-factor_refine \
		 bin/t-femtocas \
		 bin/t-expressions \
		 bin/t-generators

CHECKS=check_factor_refine \
       check_femtocas \
       check_expressions \
       check_generators

all: bin/arbtkf91

# force the test scripts to be built by listing them as requirements
tests: $(TEST_EXECUTABLES)

# run the test scripts
check: $(CHECKS)


# factor refine

check_factor_refine: bin/t-factor_refine
	$(VALGRIND) bin/t-factor_refine

factor_refine.o: factor_refine.c
	$(CC) factor_refine.c -c $(CFLAGS) -lflint -lgmp

bin/t-factor_refine: tests/t-factor_refine.c factor_refine.o
	$(CC) tests/t-factor_refine.c factor_refine.o \
		-o bin/t-factor_refine \
		-I. $(CFLAGS) -lflint -lgmp


# femtocas

check_femtocas: bin/t-femtocas
	$(ARB_LD_LIBRARY) $(VALGRIND) valgrind bin/t-femtocas

femtocas.o: femtocas.c
	$(CC) femtocas.c -c $(ARB_INCLUDES) $(CFLAGS) -lflint -lgmp -larb

bin/t-femtocas: tests/t-femtocas.c femtocas.o
	$(CC) tests/t-femtocas.c femtocas.o \
		-o bin/t-femtocas \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# expressions

check_expressions: bin/t-expressions
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/t-expressions

expressions.o: expressions.c
	$(CC) expressions.c -c $(ARB_INCLUDES) $(CFLAGS) femtocas.o \
		-lflint -lgmp -larb

bin/t-expressions: tests/t-expressions.c expressions.o \
	femtocas.o tkf91_rationals.o
	$(CC) tests/t-expressions.c expressions.o \
		femtocas.o tkf91_rationals.o \
		-o bin/t-expressions \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# generators

check_generators: bin/t-generators
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/t-generators

generators.o: generators.c
	$(CC) generators.c -c $(ARB_INCLUDES) $(CFLAGS) \
		femtocas.o expressions.o \
		-lflint -lgmp -larb

bin/t-generators: tests/t-generators.c generators.o \
	tkf91_generators.o tkf91_rationals.o
	$(CC) tests/t-generators.c \
		generators.o femtocas.o expressions.o \
		tkf91_generators.o tkf91_rationals.o \
		-o bin/t-generators \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# more generators stuff

rgenerators.o: rgenerators.c \
	femtocas.o expressions.o generators.o
	$(CC) rgenerators.c -c $(ARB_INCLUDES) $(CFLAGS) \
		femtocas.o expressions.o generators.o \
		-lflint -lgmp -larb

tkf91_rgenerators.o: tkf91_rgenerators.c \
	femtocas.o expressions.o generators.o rgenerators.o
	$(CC) tkf91_rgenerators.c -c $(ARB_INCLUDES) $(CFLAGS) \
		femtocas.o expressions.o rgenerators.o \
		-lflint -lgmp -larb

tkf91_generators.o: tkf91_generators.c \
	femtocas.o expressions.o generators.o
	$(CC) tkf91_generators.c -c $(ARB_INCLUDES) $(CFLAGS) \
		femtocas.o expressions.o generators.o \
		-lflint -lgmp -larb


# double precision wavefront implementation

wavefront_double.o: wavefront_double.c
	$(CC) wavefront_double.c -c $(CFLAGS) \
		-lflint -lgmp


# integer vector wavefront implementation

wavefront_hermite.o: wavefront_hermite.c
	$(CC) wavefront_hermite.c -c $(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb


# the matrix of breadcrumbs for the traceback stage of dynamic programming

breadcrumbs.o: breadcrumbs.c
	$(CC) breadcrumbs.c -c $(CFLAGS) \
		-lflint -lgmp


# arbitrary and double precision tkf91 dynamic programming

tkf91_generator_vecs.o: tkf91_generator_vecs.c
	$(CC) tkf91_generator_vecs.c \
		-c $(CFLAGS) \
		-lflint -lgmp

tkf91_dp_d.o: tkf91_dp_d.c \
	breadcrumbs.o
	$(CC) tkf91_dp_d.c \
		breadcrumbs.o \
		-c $(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb

tkf91_dp_r.o: tkf91_dp_r.c \
	breadcrumbs.o tkf91_generator_vecs.o
	$(CC) tkf91_dp_r.c \
		breadcrumbs.o tkf91_generator_vecs.o \
		-c $(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb

tkf91_dp_bound.o: tkf91_dp_bound.c \
	breadcrumbs.o
	$(CC) tkf91_dp_bound.c \
		breadcrumbs.o \
		-c $(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb


# the main executable

bin/arbtkf91: arbtkf91.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_double.o wavefront_hermite.o breadcrumbs.o \
	tkf91_generator_vecs.o \
	tkf91_dp_d.o tkf91_dp_r.o tkf91_dp_bound.o
	$(CC) arbtkf91.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_double.o wavefront_hermite.o breadcrumbs.o \
		tkf91_generator_vecs.o \
		tkf91_dp_d.o tkf91_dp_r.o tkf91_dp_bound.o \
		-o bin/arbtkf91 \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb -lm


# run an example

example1a: bin/arbtkf91
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/arbtkf91 \
		--precision arbitrary \
		--trace 1 \
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
		--precision bound \
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
		--precision bound \
		--trace 1 \
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
		--precision bound \
		--sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA \
		--sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA \
		--lambda-num 1 --lambda-den 1 \
		--mu-num 2 --mu-den 1 \
		--tau-num 1 --tau-den 10 \
		--pa-num 27 --pa-den 100 \
		--pc-num 24 --pc-den 100 \
		--pg-num 26 --pg-den 100 \
		--pt-num 23 --pt-den 100

example2c: bin/arbtkf91
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/arbtkf91 \
		--precision bound \
		--sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA \
		--sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA \
		--lambda-num 1 --lambda-den 1 \
		--mu-num 2 --mu-den 1 \
		--tau-num 1 --tau-den 10 \
		--pa-num 1 --pa-den 2 \
		--pc-num 1 --pc-den 4 \
		--pg-num 1 --pg-den 8 \
		--pt-num 1 --pt-den 8

example2d: bin/arbtkf91
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/arbtkf91 \
		--precision bound \
		--sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA \
		--sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA \
		--lambda-num 1 --lambda-den 1 \
		--mu-num 3 --mu-den 1 \
		--tau-num 1 --tau-den 10 \
		--pa-num 27 --pa-den 100 \
		--pc-num 24 --pc-den 100 \
		--pg-num 26 --pg-den 100 \
		--pt-num 23 --pt-den 100



clean:
	rm -f *.o
	rm -f bin/arbtkf91
	rm -f bin/t-factor_refine
	rm -f bin/t-femtocas
	rm -f bin/t-expressions
	rm -f bin/t-generators
