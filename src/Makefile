CC = mpic++


#DEFINE flags
# -DUSE_MPI  turns on MPI parallelization
#  DEFINE = #-DUSE_MPI -DNO_C11

CFLAGS =-O3 -g -DUSE_MPI  

LAMMPS_PATH = 
PROPhet_PATH = $(CURDIR)
LDFLAGS = -lm

SOURCEDIR=




ARCHIVE = ar
ARCHFLAGS = -rc
LAMMPS_LIB = libpair_PROPhet.a
ALL_LIB = libPROPhet.a



SOURCE = Structure.cpp Error.cpp  Tables.cpp Functional.cpp Potential.cpp \
	System_data.cpp Setup.cpp Parallel.cpp Grid_data.cpp System.cpp VASP.cpp\
	QE.cpp FHIAIMS.cpp Optimizer.cpp Functional_params.cpp xml_reader.cpp \
	pugixml.cpp NN.cpp NN_node.cpp NF.cpp NF_node.cpp Atom.cpp Analysis.cpp



OBJECTS = $(SOURCE:%.cpp=%.o)



all: $(OBJECTS)
	$(MAKE) PROPhet
	$(MAKE) libPROPhet.a
ifdef LAMMPS_PATH
	$(MAKE) libPROPhet_lammps.a
endif


# Just build PROPhet
PROPhet: $(OBJECTS)
	$(CC) $(DEFINE) $(CFLAGS) $(OBJECTS) main.cpp -o PROPhet $(LDFLAGS) 

# Build PROPhet as an external static library
libPROPhet.a: PROPhet
	$(ARCHIVE) $(ARCHFLAGS) $(ALL_LIB) $(OBJECTS)

# Build the PROPhet library for linking with the LAMMPS MD code
libPROPhet_lammps.a: $(OBJECTS) 
	$(CC) $(DEFINE) $(CFLAGS) -UNOT_MD -UUSE_MPI -c Parallel.cpp
	$(CC) $(DEFINE) $(CFLAGS) -UNOT_MD -UUSE_MPI -I$(LAMMPS_PATH) -c pair_nn.cpp -o pair_nn.o
	$(ARCHIVE) $(ARCHFLAGS) libPROPhet_lammps.a pair_nn.o $(OBJECTS)
	@touch Parallel.cpp


# Install everything
.PHONY: install
install:
ifneq ($(wildcard PROPhet),)
	cp PROPhet ../bin/
endif
ifneq ($(wildcard libPROPhet_lammps.a),)
# First, copy the necessary header file to the lammps directory
# Error out if we can't write there
	cp pair_nn.h $(LAMMPS_PATH) || ($(MAKE) lammps_install_error);
# Now, edit the LAMMPS_INCLUDE/Makefile.package file to tell
# it where to find the PROPhet library.
#       Make a backup of the Makefile.packag file or error out
	cp $(LAMMPS_PATH)/Makefile.package $(LAMMPS_PATH)/Makefile.package.backup || ($(MAKE) lammps_install_error);
	cat $(LAMMPS_PATH)/Makefile.package.backup | \
	sed 's|-I$(PROPhet_PATH)||'|sed '/PKG_INC/s|$$| -I$(PROPhet_PATH)|'  | \
	sed 's|-L$(PROPhet_PATH)||'|sed '/PKG_PATH/s|$$| -L$(PROPhet_PATH)|' | \
	sed 's|-lPROPhet_lammps||'|sed '/PKG_LIB/s|$$| -lPROPhet_lammps|'  \
	> $(LAMMPS_PATH)/Makefile.package || ($(MAKE) lammps_install_error);
	@echo;echo
	@echo "PROPhet successfully installed in LAMMPS package."
	@echo "You should now (re)-compile LAMMPS to link"
	@echo;echo
endif


.cpp.o:
	$(CC) $(DEFINE) $(INCLUDE) $(CFLAGS) -DNOT_MD -c $< -o $@

.PHONY: lammps_install_error
lammps_install_error:
	@echo;
	@echo "Could not automatically install PROPhet in LAMMPS"
	@echo "Try using \"sudo make install\" or manually link"
	@echo "To manually link PROPhet with LAMMPS you must:"
	@echo "1. copy pair_nn.h to $(LAMMPS_PATH)"
	@echo "2. In $(LAMMPS_PATH)/Makefile.package:"
	@echo "    add -I$(PROPhet_PATH)  to PKG_INC"
	@echo "    add -L$(PROPhet_PATH)  to PKG_PATH"
	@echo "    add -lPROPhet_lammps  to PKG_LIB"
	@echo
	@echo "See documentation for more details";echo;
	exit 1;


.PHONY: clean
clean:
	rm -f *.o PROPhet libpair_NN.a libPROPhet.a


