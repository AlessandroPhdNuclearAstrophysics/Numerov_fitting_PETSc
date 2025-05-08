-include /home/alessandro/petsc/petscdir.mk

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

.PHONY = all clean

SRCDIR := ../src
BINDIR := ../bin

NAME := fit_potential

TARGET := $(BINDIR)/${NAME}.x 

SOURCES_F90 := 	\
								$(SRCDIR)/libs/scattering_single_channel.f90 \
						    $(SRCDIR)/libs/numerov.f90 \
						    $(SRCDIR)/libs/spherical_bessel.f90 \
						    $(SRCDIR)/libs/linear_and_quadratic_regression.f90 \
						    $(SRCDIR)/libs/optimization.f90 \
						    $(SRCDIR)/libs/potential_pw.f90 \
						    $(SRCDIR)/libs/potentials/EFT_pless_new.f90 \
						    $(SRCDIR)/libs/potentials/EFT_pless_fit.f90 \

								
SOURCES_F   :=  $(SRCDIR)/libs/potentials/av14.f \
								$(SRCDIR)/libs/potentials/av18.f \
								$(SRCDIR)/libs/potentials/EFT_pless.f \
								$(SRCDIR)/libs/potentials/EFT_pless_utils.f

SOURCES_C 	:=  $(SRCDIR)/${NAME}.c
HEADERS			:=  $(SRCDIR)/libs/scattering_single_channel.h \
								$(SRCDIR)/libs/physical_constants.h
						

OBJECTS_F90 := $(patsubst $(SRCDIR)/%.f90, %.o, $(SOURCES_F90))
OBJECTS_F   := $(patsubst $(SRCDIR)/%.f, %.o, $(SOURCES_F))
OBJECTS_C   := $(patsubst $(SRCDIR)/%.c, %.o, $(SOURCES_C))

RM = rm -vf

LIBS += -L/cern/2006b/lib -llapack -lblas -lgsl -lgslcblas -fno-automatic -finit-local-zero 
# FFLAGS += -w -fno-automatic -fdefault-double-8 -fdefault-real-8 -O4
FFLAGS += -w -fno-automatic -fdefault-double-8 -fdefault-real-8 -g \
				 -ffpe-trap=denormal,invalid,zero,overflow -fcheck=all -finit-real=snan -ffpe-trap=invalid

all:	$(TARGET)
 
$(OBJECTS_F90): %.o: $(SRCDIR)/%.f90 
	${FC} $< $(FFLAGS) -c -o $@

$(OBJECTS_F): %.o: $(SRCDIR)/%.f 
	${FC} $< $(FFLAGS) -c -o $@

$(OBJECTS_C): %.o: $(SRCDIR)/%.c $(HEADERS)
	${PETSC_COMPILE_SINGLE} $< $(CFLAGS) -c -o $@
	
$(TARGET): $(OBJECTS_F90) $(OBJECTS_F) $(OBJECTS_C)
	${CLINKER} $^ $(FFLAGS) $(LIBS) $(LDLIBS) -o $@

CLEANFILES += $(OBJECTS_F90) $(OBJECTS_F) $(OBJECTS_C) $(TARGET)
