## Process this file with automake to produce Makefile.in

## For Makefile.am noobs like me, see
## https://www.gnu.org/software/automake/manual/html_node/General-Operation.html
## https://www.lrde.epita.fr/~adl/dl/autotools.pdf


bin_PROGRAMS = arbtkf91-align arbtkf91-check arbtkf91-image arbtkf91-bench
check_PROGRAMS = t-expressions t-factor_refine t-femtocas t-generators
TESTS = $(check_PROGRAMS)

CORE_SOURCES =  \
	bound_mat.c \
	count_solutions.c \
	expressions.c \
	factor_refine.c \
	femtocas.c \
	generators.c \
	model_params.c \
	rgenerators.c \
	tkf91_dp_bound.c \
	tkf91_dp.c \
	tkf91_dp_d.c \
	tkf91_dp_f.c \
	tkf91_dp_r.c \
	tkf91_generators.c \
	tkf91_generator_vecs.c \
	tkf91_rationals.c \
	tkf91_rgenerators.c \
	vis.c \
	wavefront_hermite.c \
	bound_mat.h \
	count_solutions.h \
	expressions.h \
	factor_refine.h \
	femtocas.h \
	generators.h \
	model_params.h \
	printutil.h \
	rgenerators.h \
	tkf91_dp_bound.h \
	tkf91_dp_d.h \
	tkf91_dp_f.h \
	tkf91_dp.h \
	tkf91_dp_r.h \
	tkf91_generator_indices.h \
	tkf91_generators.h \
	tkf91_generator_vecs.h \
	tkf91_rationals.h \
	tkf91_rgenerators.h \
	vis.h \
	wavefront_hermite.h \
	dp.c \
	dp.h \
	forward.c \
	forward.h \
	unused.h

JSON_SOURCES = \
	json_model_params.c \
	jsonutil.c \
	runjson.c \
	json_model_params.h \
	jsonutil.h \
	runjson.h

ALL_SOURCES = $(CORE_SOURCES) $(JSON_SOURCES)

t_expressions_SOURCES =  $(CORE_SOURCES) t-expressions.c
t_factor_refine_SOURCES =  $(CORE_SOURCES) t-factor_refine.c
t_femtocas_SOURCES =  $(CORE_SOURCES) t-femtocas.c
t_generators_SOURCES =  $(CORE_SOURCES) t-generators.c

arbtkf91_align_SOURCES =  $(ALL_SOURCES) arbtkf91-align.c
arbtkf91_bench_SOURCES =  $(ALL_SOURCES) arbtkf91-bench.c
arbtkf91_image_SOURCES =  $(ALL_SOURCES) arbtkf91-image.c
arbtkf91_check_SOURCES =  $(ALL_SOURCES) arbtkf91-check.c
