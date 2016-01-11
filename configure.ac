AC_INIT([amarbtkf91], [1.0], [bug-report@address])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])

## https://autotools.io/autoconf/arguments.html
## https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html

# But also see http://stackoverflow.com/questions/10199754



AC_ARG_WITH([flint],
	    [AS_HELP_STRING([--with-flint],
			    [Specify location of flint @<:@default=check@:>@])],
	    [],
	    [with_flint=check])

AC_ARG_WITH([arb],
	    [AS_HELP_STRING([--with-arb],
			    [Specify location of arb @<:@default=check@:>@])],
	    [],
	    [with_arb=check])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
