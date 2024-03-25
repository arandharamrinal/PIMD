FCOMP= ifort
OPT= -O4 
#For compilattion without PLUMED 
FFLAGS= -r8 -mavx2 -mkl -qopenmp #-fpp  $(PLUMED_LOAD) 
#For compilattion with PLUMED 
#FFLAGS= -r8 -mavx2 -mkl -qopenmp -I$(VSLDIR) $(PLUMED_LOAD)  
#Intel mkl library 
SRCDIR=./src
BINDIR=./bin
BUILDDIR=./build
VSLDIR=/home/mrinal/intel/compilers_and_libraries_2019.5.281/linux/mkl/include
VSLFILE=$(VSLDIR)/mkl_vsl.f90
FILES=  \
    mkl_vsl.f90\
    UnitConversion.f90\
    allvariables.f90\
	brngvars.f90\
    USvariables.f90\
	gaussModule.f90\
    plumedModule.f90\
	sgdmltofortran.f90\
    house.f90\
    allocatearrays.f90\
    createfiles.f90\
    checkVariables.f90\
    initializeVariables.f90\
    nmmassmat.f90\
    rpolymat.f90\
    nmmatrix.f90\
    nmtrans.f90\
    genrn.f90\
    thermostats_vars.f90\
    umbrellaSampling2D.f90\
    getGaussForce.f90\
    computeProperties.f90\
    computeSysEnergy.f90\
    computeThermostatEnergy.f90\
    forceVelCorrection.f90\
    integrateMNHC.f90\
    integratePIGLET.f90\
    LangevinDynamics.f90\
    writeTofiles.f90\
    freeRingPolymerUpdate.f90\
    integrateClassicalMD.f90\
    integrateClassicalNVE.f90\
    integrateRPMD.f90\
    integrateTRPMD.f90\
    integratePIMDPC.f90\
    integratePIMDNPC.f90\
    setupClassicalNVE.f90\
    setupClassicalMD.f90\
    setupTRPMD.f90\
    setupPIMDPC.f90\
    setupPIMDNPC.f90\
    pathIntegral.f90\
#    plumedInitialize.f90\
#    plumedForceCalc.f90\

FFILES=$(addprefix $(SRCDIR)/,$(FILES))
OBJFILES=$(patsubst %.f90,$(BUILDDIR)/%.o,$(FILES))
BINARY=$(BINDIR)/pimd.exe

.PHONY: all clean

all: $(BINARY)

$(BINARY): $(OBJFILES)
	mkdir -p $(BINDIR)
	$(FCOMP) $(FFLAGS) $^ -o $@ 

$(BUILDDIR)/%.o:$(SRCDIR)/%.f90
	mkdir -p $(BUILDDIR)
	$(FCOMP) $(OPT) $(FFLAGS) -c $^ -o $@ -module $(BUILDDIR)

clean:
	rm -rf $(BINDIR) $(OBJFILES) $(BINARY)

