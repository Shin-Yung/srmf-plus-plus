#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
#include ../root_dir.mk # must be included first
include ./make_options.mk
#-------------------------------------------------------------
# Source files
SRCS = scheduler/mpi_comm.cpp 
SRCS+= scheduler/cmdargs.cpp 
SRCS+= scheduler/inputparams.cpp 
SRCS+= scheduler/taskparams.cpp 
SRCS+= scheduler/worker.cpp 
SRCS+= scheduler/master_scheduler.cpp
SRCS+= scheduler/scheduler.cpp
#SRCS+= xml/pugixml.cpp 
#SRCS+= expression/expression.cpp 
#SRCS+= expression/tokens.cpp 
#SRCS+= expression/functions.cpp
#SRCS+= expression/objects.cpp
SRCS+= expression/complex_expression.cpp
SRCS+= utils/utils.cpp 
SRCS+= lattice/lattice.cpp
SRCS+= lattice/latticelibrary.cpp
SRCS+= lattice/graph.cpp
SRCS+= lattice/blochbasis.cpp
SRCS+= model/strmatrix.cpp
SRCS+= model/hamiltonian_term.cpp
SRCS+= model/model.cpp
SRCS+= model/modellibrary.cpp
SRCS+= srmf/spinon.cpp
SRCS+= srmf/srmf.cpp
SRCS+= main.cpp
VMC_SRCS = $(addprefix src/,$(SRCS))
#-------------------------------------------------------------
# Headers
HDRS=    scheduler/mpi_comm.h \
         scheduler/optionparser.h scheduler/cmdargs.h \
         scheduler/inputparams.h scheduler/worker.h scheduler/task.h \
         scheduler/scheduler.h \
         utils/complex_expression.h \
         utils/utils.h \
         lattice/constants.h lattice/lattice.h lattice/graph.h \
	 lattice/blochbasis.h \
	 montecarlo/simulator.h \
         model/strmatrix.h model/modelparams.h  model/quantum_op.h \
	 model/hamiltonian_term.h \
	 srmf/spinon.h \
	 model/model.h \
	 srmf/srmf.h 
#         expression/expression.h expression/shunting_yard.h \
         expression/tokens.h expression/functions.h expression/objects.h \
         expression/pack.h \
VMC_HDRS = $(addprefix src/,$(HDRS))
MUPARSER_LIB = $(PROJECT_ROOT)/src/expression/muparserx/libmuparserx.a
#-------------------------------------------------------------
# Target
TAGT=a.out
ifeq ($(MPI), HAVE_BOOST_MPI)
TAGT=v.out
endif

# Put all auto generated stuff to this build dir.
ifeq ($(BUILD_DIR), $(CURDIR))
  $(error In-source build is not allowed, choose another build directory)
endif

# All .o files go to BULD_DIR
OBJS=$(patsubst %.cpp,$(BUILD_DIR)/%.o,$(VMC_SRCS))
# GCC/Clang will create these .d files containing dependencies.
DEPS=$(patsubst %.o,%.d,$(OBJS)) 
# compiler flags

.PHONY: all
all: $(TAGT) #$(INCL_HDRS)

$(TAGT): $(OBJS)
	$(VMC_CXX) -o $(TAGT) $(OBJS) $(VMC_LDFLAGS) $(VMC_LIBS) $(MUPARSER_LIB) 

%.o: %.cpp
	$(VMC_CXX) -c $(VMC_CXXFLAGS) -o $@ $<

# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	@echo "$(VMC_CXX) -c $(VMC_CXXFLAGS) -o $(@F) $(<F)"
	@$(VMC_CXX) -MMD -c $(VMC_CXXFLAGS) -o $@ $<

$(VMC_INCLDIR)/%.h: %.h 
	@mkdir -p $(@D)
	@echo "Copying $< to 'include'" 
	@cp -f $< $@

# installation
#prefix = ../install#/usr/local
#libdir = $(prefix)/lib
#includedir = $(prefix)/include/cmc++

.PHONY: install
install:	
	@echo "Already installed in $(VMC_LIBDIR) and $(VMC_INCLDIR)" 

.PHONY: clean
clean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
	@echo "Removing $(TAGT)"
	@rm -f $(TAGT) 

.PHONY: bclean
bclean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
