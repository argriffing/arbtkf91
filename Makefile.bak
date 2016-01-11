
# work around https://github.com/fredrik-johansson/arb/issues/24
ARB_INCLUDES=-I/usr/local/include/flint

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

MAIN_EXECUTABLES=bin/arbtkf91-check \
		bin/arbtkf91-bench \
		bin/arbtkf91-align \
		bin/arbtkf91-image

TEST_EXECUTABLES=bin/t-factor_refine \
		 bin/t-femtocas \
		 bin/t-expressions \
		 bin/t-generators

CHECKS=check_factor_refine \
       check_femtocas \
       check_expressions \
       check_generators \
       python_tests

all: $(MAIN_EXECUTABLES)

# force the test scripts to be built by listing them as requirements
tests: $(TEST_EXECUTABLES)

# run the test scripts
check: $(CHECKS)


# run the python tests

python_tests: $(MAIN_EXECUTABLES)
	nosetests

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
	$(VALGRIND) valgrind bin/t-femtocas

femtocas.o: femtocas.c
	$(CC) femtocas.c -c $(ARB_INCLUDES) $(CFLAGS)

bin/t-femtocas: tests/t-femtocas.c femtocas.o
	$(CC) tests/t-femtocas.c femtocas.o \
		-o bin/t-femtocas \
		-I. $(ARB_INCLUDES) $(CFLAGS) -lflint -lgmp -larb


# expressions

check_expressions: bin/t-expressions
	$(VALGRIND) bin/t-expressions

expressions.o: expressions.c \
	femtocas.h
	$(CC) expressions.c -c $(ARB_INCLUDES) $(CFLAGS)

bin/t-expressions: tests/t-expressions.c expressions.o \
	femtocas.o tkf91_rationals.o
	$(CC) tests/t-expressions.c expressions.o \
		femtocas.o tkf91_rationals.o \
		-o bin/t-expressions \
		-I. $(ARB_INCLUDES) $(CFLAGS) -lflint -lgmp -larb


# generators

check_generators: bin/t-generators
	$(VALGRIND) bin/t-generators

generators.o: generators.c \
	femtocas.h expressions.h
	$(CC) generators.c -c $(ARB_INCLUDES) $(CFLAGS)

bin/t-generators: tests/t-generators.c generators.o \
	tkf91_generators.o tkf91_rationals.o
	$(CC) tests/t-generators.c \
		generators.o femtocas.o expressions.o \
		tkf91_generators.o tkf91_rationals.o \
		-o bin/t-generators \
		-I. $(ARB_INCLUDES) $(CFLAGS) -lflint -lgmp -larb


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
	$(CC) wavefront_hermite.c -c $(ARB_INCLUDES) $(CFLAGS)


# the matrix of breadcrumbs for the traceback stage of dynamic programming

breadcrumbs.o: breadcrumbs.c breadcrumbs.h
	$(CC) breadcrumbs.c -c $(CFLAGS)

vis.o: vis.c breadcrumbs.h
	$(CC) vis.c -c $(CFLAGS)

# arbitrary and double precision tkf91 dynamic programming

tkf91_generator_vecs.o: tkf91_generator_vecs.c
	$(CC) tkf91_generator_vecs.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp_d.o: tkf91_dp_d.c tkf91_dp_d.h \
	breadcrumbs.h printutil.h tkf91_dp.h
	$(CC) tkf91_dp_d.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp_f.o: tkf91_dp_f.c tkf91_dp_f.h \
	printutil.h tkf91_dp.h
	$(CC) tkf91_dp_f.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp_r.o: tkf91_dp_r.c tkf91_dp_r.h \
	breadcrumbs.h tkf91_generator_vecs.h tkf91_dp.h
	$(CC) tkf91_dp_r.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp_bound.o: tkf91_dp_bound.c \
	breadcrumbs.h bound_mat.h printutil.h vis.h count_solutions.h tkf91_dp.h
	$(CC) tkf91_dp_bound.c \
		-c $(ARB_INCLUDES) $(CFLAGS)

bound_mat.o: bound_mat.c breadcrumbs.h
	$(CC) bound_mat.c -c $(ARB_INCLUDES) $(CFLAGS)

tkf91_dp.o: tkf91_dp.c tkf91_dp.h
	$(CC) tkf91_dp.c -c $(ARB_INCLUDES) $(CFLAGS)

count_solutions.o: count_solutions.c count_solutions.h breadcrumbs.h
	$(CC) count_solutions.c -c $(CFLAGS)


# json related stuff

runjson.o: runjson.c runjson.h
	$(CC) runjson.c -c $(CFLAGS)

jsonutil.o: jsonutil.c jsonutil.h
	$(CC) jsonutil.c -c $(CFLAGS)

model_params.o: model_params.c model_params.h
	$(CC) model_params.c -c $(CFLAGS)

json_model_params.o: json_model_params.c json_model_params.h
	$(CC) json_model_params.c -c $(CFLAGS)



# the main executables

bin/arbtkf91-check: arbtkf91-check.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_hermite.o breadcrumbs.o \
	tkf91_dp.o tkf91_generator_vecs.o count_solutions.o \
	tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o bound_mat.o \
	runjson.o jsonutil.o json_model_params.o model_params.o vis.o
	$(CC) arbtkf91-check.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_hermite.o breadcrumbs.o \
		tkf91_dp.o tkf91_generator_vecs.o \
		tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o \
		bound_mat.o count_solutions.o \
		runjson.o jsonutil.o json_model_params.o model_params.o vis.o \
		-o bin/arbtkf91-check \
		$(ARB_INCLUDES) $(CFLAGS) \
		-lflint -lgmp -larb -lm -lpng -ljansson

bin/arbtkf91-bench: arbtkf91-bench.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_hermite.o breadcrumbs.o \
	tkf91_dp.o tkf91_generator_vecs.o count_solutions.o \
	tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o bound_mat.o \
	runjson.o jsonutil.o json_model_params.o model_params.o vis.o
	$(CC) arbtkf91-bench.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_hermite.o breadcrumbs.o \
		tkf91_dp.o tkf91_generator_vecs.o \
		tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o \
		bound_mat.o count_solutions.o \
		runjson.o jsonutil.o json_model_params.o model_params.o vis.o \
		-o bin/arbtkf91-bench \
		$(ARB_INCLUDES) $(CFLAGS) \
		-lflint -lgmp -larb -lm -lpng -ljansson

bin/arbtkf91-align: arbtkf91-align.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_hermite.o breadcrumbs.o \
	tkf91_dp.o tkf91_generator_vecs.o count_solutions.o \
	tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o bound_mat.o \
	runjson.o jsonutil.o json_model_params.o model_params.o vis.o
	$(CC) arbtkf91-align.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_hermite.o breadcrumbs.o \
		tkf91_dp.o tkf91_generator_vecs.o \
		tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o \
		bound_mat.o count_solutions.o \
		runjson.o jsonutil.o json_model_params.o model_params.o vis.o \
		-o bin/arbtkf91-align \
		$(ARB_INCLUDES) $(CFLAGS) \
		-lflint -lgmp -larb -lm -lpng -ljansson

bin/arbtkf91-image: arbtkf91-image.c \
	femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
	generators.o tkf91_generators.o \
	rgenerators.o tkf91_rgenerators.o \
	wavefront_hermite.o breadcrumbs.o \
	tkf91_dp.o tkf91_generator_vecs.o count_solutions.o \
	tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o bound_mat.o \
	runjson.o jsonutil.o json_model_params.o model_params.o vis.o
	$(CC) arbtkf91-image.c \
		femtocas.o factor_refine.o expressions.o tkf91_rationals.o \
		generators.o tkf91_generators.o \
		rgenerators.o tkf91_rgenerators.o \
		wavefront_hermite.o breadcrumbs.o \
		tkf91_dp.o tkf91_generator_vecs.o \
		tkf91_dp_d.o tkf91_dp_f.o tkf91_dp_r.o tkf91_dp_bound.o \
		bound_mat.o count_solutions.o \
		runjson.o jsonutil.o json_model_params.o model_params.o vis.o \
		-o bin/arbtkf91-image \
		$(ARB_INCLUDES) $(CFLAGS) \
		-lflint -lgmp -larb -lm -lpng -ljansson




clean:
	rm -f *.o
	rm -f *.pyc
	rm -f bin/arbtkf91-check
	rm -f bin/arbtkf91-bench
	rm -f bin/arbtkf91-align
	rm -f bin/arbtkf91-image
	rm -f bin/t-factor_refine
	rm -f bin/t-femtocas
	rm -f bin/t-expressions
	rm -f bin/t-generators
