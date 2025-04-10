PROGRAM= tlaloque
COMPILER= gfortran
######  USER_FLAGS for gfortran:
#USER_FLAGS= -O3 -Wall -fcheck=all
######  USER_FLAGS for ifort:
USER_FLAGS= -O3 #-warn all -check all
MPI= N
DOUBLEP= Y
PASB= N
GRV= N
COOL= N

MODULES_USER= \

MODULES_MAIN= \
./source/constants.o \
./user/parameters.o \
./source/globals.o \
./source/cooling.o \
./source/pointsource_gravity.o \
./source/sources.o \
./source/winds.o \
./user/user.o 

OBJECTS_MAIN= \
./source/HydroCore.o \
./source/init.o \
./source/boundaries.o \
./source/HLLC.o \
./source/mainHD.o 

CFLAGS = $(USER_FLAGS) -cpp

ifeq ($(MPI),Y)
CFLAGS += -DMPIP
endif

ifeq ($(COOL),Y)
CFLAGS += -DCOOL
endif

ifeq ($(GRV),Y)
CFLAGS += -DGRAV
endif

ifeq ($(DOUBLEP),Y)
CFLAGS += -DDOUBLEP
ifeq ($(COMPILER),ifort)
CFLAGS += -r8
endif
ifeq ($(COMPILER),gfortran)
CFLAGS += -fdefault-real-8
endif
endif

ifeq ($(PASB),Y)
CFLAGS += -DPASBP
endif

ifeq ($(MPI),Y)
COMPILER = mpif90
endif

OBJECTS_ALL = ${MODULES_MAIN} ${MODULES_USER} ${OBJECTS_MAIN}

$(PROGRAM) : prebuild ${OBJECTS_ALL}
	@echo Linking object files ...
	@$(COMPILER) $(CFLAGS) $(OBJECTS_ALL) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod source/*.o source/*.mod user/*.o user/*.mod
	@echo "Done! (`date`)"

prebuild :
	@echo "$(PROGRAM) build started `date`"

%.o : %.f90
	@echo Compiling $^ ...
	@$(COMPILER) $(CFLAGS) -c $^ -o $@

clean :
	rm -f *.o *.mod source/*.o source/*.mod user/*.o user/*.mod

cleanall :
	rm -f *.o *.mod source/*.o source/*.mod user/*.o user/*.mod
	rm -f $(PROGRAM).*
	rm -f coldens
	rm -f extract
	rm -f DATA/*.bin
	rm -f DATA/*.vtk
	rm -f DATA/*.dat
	rm -f DATA/*.log
	rm -f DATA/*.visit
