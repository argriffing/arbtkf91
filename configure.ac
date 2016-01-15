AC_INIT([arbtkf91], [0.0.1], [bug-report@address])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])

# https://www.gnu.org/software/autoconf-archive/ax_valgrind_check.html
#AX_VALGRIND_CHECK

# http://stackoverflow.com/questions/5298830
AC_CONFIG_MACRO_DIR([autoconf_macros])
AX_CC_MAXOPT


# Manually written autoconf code to look for library functions...

AC_SEARCH_LIBS(fmaxf, [m], [], [AC_MSG_ERROR(
	[unable to find the fmaxf() function (math library missing?)])])

AC_SEARCH_LIBS(json_equal, [jansson], [], [AC_MSG_ERROR(
	[unable to find the json_equal() function (jansson library missing?)])])

AC_SEARCH_LIBS(png_create_write_struct, [png], [], [AC_MSG_ERROR(
	[unable to find the png_create_write_struct() function (png library missing?)])])

AC_SEARCH_LIBS(__gmpz_set, [gmp], [], [AC_MSG_ERROR(
	[unable to find the __gmpz_set() function (gmp library missing?)])])

AC_SEARCH_LIBS(fmpz_init, [flint], [], [AC_MSG_ERROR(
	[unable to find the fmpz_init() function (flint library missing?)])])

AC_SEARCH_LIBS(arb_init, [arb], [], [AC_MSG_ERROR(
	[unable to find the arb_init() function (arb library missing?)])])


AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
