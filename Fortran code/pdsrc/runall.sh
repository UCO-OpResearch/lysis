#!/bin/bash
################################################################################
GCC(){
exec 2>&1
set -x
gcc -c macros.c  # build auxiliary C routines
gfortran -J. -c ncurses.f90
gfortran -c nc_printhtml.f90 nc_errmessage.f90
for NAME in basics/*.f90
do
   gfortran -I. -J. $NAME nc_printhtml.o nc_errmessage.o ncurses.o macros.o -lncurses -o $(basename $NAME .f90)
   ./$(basename $NAME .f90)
done
}
################################################################################
INTEL(){
exec 2>&1
set -x
icc -c macros.c  # build auxiliary C routines
ifort -I. -c ncurses.f90
ifort -I. -c nc_printhtml.f90 nc_errmessage.f90
for NAME in basics/*.f90
do
   ifort -I. $NAME nc_printhtml.o nc_errmessage.o ncurses.o macros.o -lncurses -o $(basename $NAME .f90)
   ./$(basename $NAME .f90)
done
}
################################################################################
if [ "$(which ifort 2>/dev/null)" != ' ' ]
then
   INTEL
else
   GCC
fi
################################################################################
exit
################################################################################
