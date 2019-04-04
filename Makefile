TARGET = a.out 
OBJECTS = toukei.o MTlt1.o boot.o global.o \
	inisetting.o syscalcSG.o metropolis.o hcSGXY.o
MOD_FILES  = global.mod iniset.mod metropolis_method.mod\
	mtmod.mod syscalc.mod toukei.mod boot.mod
#FC = gfortran
FC = ifort
#FC = mpif90

FFLAGS =
LDFLAGS =


# for gfortran
ifeq (${FC},gfortran)
#	FFLAGS += -fimplicit-none 
	LDFLAGS += -fopenmp -fbounds-check 
endif

# for ifort
ifeq (${FC},ifort)
#	MKLROOT = /opt/intel/composer_xe_2011_sp1.7.256/mkl
#	FFLAGS += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include
#	LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmkl_blas95.a
#	LDFLAGS += ${MKLROOT}/lib/ia32/libmkl_lapack95.a
	LDFLAGS +=  -O3 -openmp
endif

# for mpif90
ifeq (${FC},mpif90)
#       MKLROOT = /opt/intel/composer_xe_2011_sp1.7.256/mkl
#       FFLAGS += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include
#       LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmkl_blas95.a
#       LDFLAGS += ${MKLROOT}/lib/ia32/libmkl_lapack95.a
        LDFLAGS += -qopenmp -O3 -ipo
endif



.SUFFIXES : .o .f90
.f90.o:
	${FC} -c ${LDFLAGS} $<

${TARGET} : ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS} 


.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
