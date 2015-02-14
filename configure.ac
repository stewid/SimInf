# siminf, a framework for stochastic disease spread simulations
# Copyright (C) 2015  Pavol Bauer
# Copyright (C) 2015  Stefan Engblom
# Copyright (C) 2015  Stefan Widgren
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.69])
AC_INIT([siminf], [see.DESCRIPTION.file], [https://github.com/stewid/siminf/issues])
AC_CONFIG_SRCDIR([src/siminf.c])
AC_CONFIG_AUX_DIR([tools])

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC

# Find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

LIBS="-lm -lgsl -lgslcblas ${LIBS}"

AC_SUBST(LIBS)

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
