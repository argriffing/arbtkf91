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


CFLAGS+= $(WARNINGS) -O3 $(DEFS) $(GLOBAL_DEFS) -march=native -ffast-math
CFLAGS+= -g
#CFLAGS+= -std=c++11


TEST_EXECUTABLES=bin/t-factor_refine \
		 bin/t-femtocas \
		 bin/t-expressions \
		 bin/t-generators

CHECKS=check_factor_refine \
       check_femtocas \
       check_expressions \
       check_generators \
       check_arbtkf91_check

all: bin/arbtkf91 bin/arbtkf91-check bin/arbtkf91-bench

# force the test scripts to be built by listing them as requirements
tests: $(TEST_EXECUTABLES)

# run the test scripts
check: $(CHECKS)


# a regression test

check_arbtkf91_check: bin/arbtkf91-check
	bin/arbtkf91-check < tests/data/a_in.json > tmp.json
	python tests/compare_json.py tests/data/a_out.json tmp.json
	#
	bin/arbtkf91-check < tests/data/b_in.json > tmp.json
	python tests/compare_json.py tests/data/b_out.json tmp.json

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
	$(CC) femtocas.c -c $(ARB_INCLUDES) $(CFLAGS)

bin/t-femtocas: tests/t-femtocas.c femtocas.o
	$(CC) tests/t-femtocas.c femtocas.o \
		-o bin/t-femtocas \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# expressions

check_expressions: bin/t-expressions
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/t-expressions

expressions.o: expressions.c \
	femtocas.h
	$(CC) expressions.c -c $(ARB_INCLUDES) $(CFLAGS)

bin/t-expressions: tests/t-expressions.c expressions.o \
	femtocas.o tkf91_rationals.o
	$(CC) tests/t-expressions.c expressions.o \
		femtocas.o tkf91_rationals.o \
		-o bin/t-expressions \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# generators

check_generators: bin/t-generators
	$(ARB_LD_LIBRARY) $(VALGRIND) bin/t-generators

generators.o: generators.c \
	femtocas.h expressions.h
	$(CC) generators.c -c $(ARB_INCLUDES) $(CFLAGS)

bin/t-generators: tests/t-generators.c generators.o \
	tkf91_generators.o tkf91_rationals.o
	$(CC) tests/t-generators.c \
		generators.o femtocas.o expressions.o \
		tkf91_generators.o tkf91_rationals.o \
		-o bin/t-generators \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) -lflint -lgmp -larb


# more generators stuff

rgenerators.o: rgenerators.c \
	femtocas.h expressions.h generators.h
	$(CC) rgenerators.c -c $(ARB_INCLUDES) $(CFLAGS)

tkf91_rgenerators.o: tkf91_rgenerators.c \
	femtocas.h expressions.h generators.h rgenerators.h
	$(CC) tkf91_rgenerators.c -c $(ARB_INCLUDES) $(CFLAGS)

tkf91_generators.o: tkf91_generators.c \
	femtocas.h expressions.h generators.h
	$(CC) tkf91_generators.c -c $(ARB_INCLUDES) $(CFLAGS)


# integer vector wavefront implementation

wavefront_hermite.o: wavefront_hermite.c
	$(CC) wavefront_hermite.c -c $(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS)


# the matrix of breadcrumbs for the traceback stage of dynamic programming

breadcrumbs.o: breadcrumbs.c
	$(CC) breadcrumbs.c -c $(CFLAGS)

vis.o: vis.c breadcrumbs.h
	$(CC) vis.c -c $(CFLAGS)

# arbitrary and double precision tkf91 dynamic programming

tkf91_generator_vecs.o: tkf91_generator_vecs.c
	$(CC) tkf91_generator_vecs.c \
		-c $(CFLAGS)

tkf91_dp_d.o: tkf91_dp_d.c \
	breadcrumbs.h printutil.h tkf91_dp.h
	$(CC) tkf91_dp_d.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp_f.o: tkf91_dp_f.c \
	printutil.h tkf91_dp.h
	$(CC) tkf91_dp_f.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp_r.o: tkf91_dp_r.c \
	breadcrumbs.h tkf91_generator_vecs.h tkf91_dp.h
	$(CC) tkf91_dp_r.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp_bound.o: tkf91_dp_bound.c \
	breadcrumbs.h bound_mat.h printutil.h vis.h count_solutions.h tkf91_dp.h
	$(CC) tkf91_dp_bound.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

bound_mat.o: bound_mat.c breadcrumbs.h
	$(CC) bound_mat.c -c $(CFLAGS)

tkf91_dp.o: tkf91_dp.c
	$(CC) tkf91_dp.c -c $(ARB_INCLUDES) $(CFLAGS)

count_solutions.o: count_solutions.c count_solutions.h breadcrumbs.h
	$(CC) count_solutions.c -c $(CFLAGS)


# a json thing

runjson.o: runjson.c runjson.h
	$(CC) runjson.c -c $(CFLAGS)

jsonutil.o: jsonutil.c jsonutil.h
	$(CC) jsonutil.c -c $(CFLAGS)



# the main executable

bin/arbtkf91: arbtkf91.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_hermite.o breadcrumbs.o \
	tkf91_dp.o tkf91_generator_vecs.o vis.o count_solutions.o \
	tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o bound_mat.o
	$(CC) arbtkf91.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_hermite.o breadcrumbs.o \
		tkf91_dp.o tkf91_generator_vecs.o vis.o \
		tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o \
		bound_mat.o count_solutions.o \
		-o bin/arbtkf91 \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb -lm -lpng

bin/arbtkf91-check: arbtkf91-check.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_hermite.o breadcrumbs.o \
	tkf91_dp.o tkf91_generator_vecs.o count_solutions.o \
	tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o bound_mat.o \
	runjson.o jsonutil.o model_params.o vis.o
	$(CC) arbtkf91-check.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_hermite.o breadcrumbs.o \
		tkf91_dp.o tkf91_generator_vecs.o \
		tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o \
		bound_mat.o count_solutions.o \
		runjson.o jsonutil.o model_params.o vis.o \
		-o bin/arbtkf91-check \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb -lm -lpng -ljansson

bin/arbtkf91-bench: arbtkf91-bench.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_hermite.o breadcrumbs.o \
	tkf91_dp.o tkf91_generator_vecs.o count_solutions.o \
	tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o bound_mat.o \
	runjson.o jsonutil.o model_params.o vis.o
	$(CC) arbtkf91-bench.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_hermite.o breadcrumbs.o \
		tkf91_dp.o tkf91_generator_vecs.o \
		tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o \
		bound_mat.o count_solutions.o \
		runjson.o jsonutil.o model_params.o vis.o \
		-o bin/arbtkf91-bench \
		$(ARB_INCLUDES) $(ARB_LIBS) $(CFLAGS) \
		-lflint -lgmp -larb -lm -lpng -ljansson


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
	rm -f bin/arbtkf91-check
	rm -f bin/arbtkf91-bench
	rm -f bin/t-factor_refine
	rm -f bin/t-femtocas
	rm -f bin/t-expressions
	rm -f bin/t-generators
