## Process this file with automake to produce Makefile.in

## For Makefile.am noobs like me, see
## https://www.gnu.org/software/automake/manual/html_node/General-Operation.html
## https://www.lrde.epita.fr/~adl/dl/autotools.pdf


bin_PROGRAMS = arbtkf91-align arbtkf91-check arbtkf91-image arbtkf91-bench

ALL_LDADD=-lflint -lgmp -larb -lm -lpng -ljansson

ALL_SOURCES =  \
	bound_mat.c \
	breadcrumbs.c \
	count_solutions.c \
	expressions.c \
	factor_refine.c \
	femtocas.c \
	generators.c \
	json_model_params.c \
	jsonutil.c \
	model_params.c \
	rgenerators.c \
	runjson.c \
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
	breadcrumbs.h \
	count_solutions.h \
	expressions.h \
	factor_refine.h \
	femtocas.h \
	generators.h \
	json_model_params.h \
	jsonutil.h \
	model_params.h \
	printutil.h \
	rgenerators.h \
	runjson.h \
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
	wavefront_hermite.h


arbtkf91_align_LDADD = $(ALL_LDADD)
arbtkf91_align_SOURCES =  $(ALL_SOURCES) arbtkf91-align.c

arbtkf91_bench_LDADD = $(ALL_LDADD)
arbtkf91_bench_SOURCES =  $(ALL_SOURCES) arbtkf91-bench.c

arbtkf91_image_LDADD = $(ALL_LDADD)
arbtkf91_image_SOURCES =  $(ALL_SOURCES) arbtkf91-image.c

arbtkf91_check_LDADD = $(ALL_LDADD)
arbtkf91_check_SOURCES =  $(ALL_SOURCES) arbtkf91-check.c
