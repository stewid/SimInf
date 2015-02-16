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

# The configuration for OpenMP is from an example by Simon Urbanek on
# the r-devel mail-list:
# https://stat.ethz.ch/pipermail/r-devel/2010-September/058488.html

AC_PREREQ([2.69])
AC_INIT([siminf], [see.DESCRIPTION.file], [https://github.com/stewid/siminf/issues])
AC_CONFIG_SRCDIR([src/siminf.c])
AC_CONFIG_AUX_DIR([tools])

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC

# Check for OpenMP
AC_OPENMP

# Find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`
CFLAGS=`"${RBIN}" CMD config CFLAGS`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
LIBS="${PKG_LIBS}"

# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp=no
if test -n "${OPENMP_CFLAGS}"; then
  AC_MSG_CHECKING([whether OpenMP will work in a package])
  AC_LANG_CONFTEST(
  [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
  PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
  AC_MSG_RESULT([${ac_pkg_openmp}])
fi

if test "x${ac_pkg_openmp}" = xyes; then
        SIMINF_HAVE_OPENMP=TRUE
        CPPFLAGS="${CPPFLAGS} -DSIMINF_OMP"
else
        OPENMP_CFLAGS=''
        SIMINF_HAVE_OPENMP=FALSE
        AC_MSG_WARN([
  -----------------------------------------------
   Unable to find OpenMP on this system.
   Building a version without support for OpenMP.
  -----------------------------------------------])
fi

# Check if SuperGlue is available
# FIXME: Add support for SuperGlue
ac_pkg_superglue=no

if test "x${ac_pkg_superglue}" = xyes; then
        SIMINF_HAVE_SUPERGLUE=TRUE
else
        SIMINF_HAVE_SUPERGLUE=FALSE
        AC_MSG_WARN([
  ---------------------------------------------------
   Unable to use SuperGlue on this system.
   Building a version without support for SuperGlue.
  ---------------------------------------------------])
fi

LIBS="-lm -lgsl -lgslcblas ${LIBS}"

AC_SUBST(SIMINF_HAVE_OPENMP)
AC_SUBST(SIMINF_HAVE_SUPERGLUE)
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)

AC_CONFIG_FILES([R/dependencies.r src/Makevars])

AC_OUTPUT
