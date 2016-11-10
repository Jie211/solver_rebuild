PROGRAM_NAME := Solver

program_CXX_SRCS := $(wildcard *.cpp)
program_CXX_SRCS += $(wildcard ./CRS/*.cpp) #Find C++ source files from additonal directories
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}

program_CU_SRCS := $(wildcard *.cu)
program_CU_OBJS := ${program_CU_SRCS:.cu=.cuo}

program_INCLUDE_DIRS := #C++ Include directories

# Compiler flags
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
# CXXFLAGS += -g -O3 -std=c++0x -Wall -pedantic
CXXFLAGS += -g -Wall -fopenmp 

GEN_SM35 := -gencode=arch=compute_35,code=\"sm_35,compute_35\" #Target CC 3.5, for example
# NVFLAGS := -O3 -rdc=true #rdc=true needed for separable compilation
NVFLAGS := -rdc=true -Xcompiler "-fopenmp"#rdc=true needed for separable compilation
NVFLAGS += $(GEN_SM35)
NVFLAGS += $(foreach includedir,$(program_CU_INCLUDE_DIRS),-I$(includedir))

CUO_O_OBJECTS := ${program_CU_OBJS:.cuo=.cuo.o}

OBJECTS = $(program_CU_OBJS) $(program_CXX_OBJS)

.PHONY: all clean distclean

all: $(PROGRAM_NAME) 

debug: CXXFLAGS = -g -O0 -std=c++0x -Wall -pedantic -DDEBUG $(EXTRA_FLAGS)
debug: NVFLAGS = -O0 $(GEN_SM35) -g -G
debug: NVFLAGS += $(foreach includedir,$(program_CU_INCLUDE_DIRS),-I$(includedir))
debug: $(PROGRAM_NAME)

%.cuo: %.cu %.cuh
	nvcc $(NVFLAGS) $(CPPFLAGS) -o $@ -dc $<

$(PROGRAM_NAME): $(OBJECTS)
	@ for cu_obj in $(program_CU_OBJS); \
	do                                  \
		mv $$cu_obj $$cu_obj.o;           \
	done
	nvcc $(NVFLAGS) $(CPPFLAGS) -o $@ $(program_CXX_OBJS) $(CUO_O_OBJECTS)
	@ for cu_obj in $(CUO_O_OBJECTS);   \
	do                                  \
		mv $$cu_obj $${cu_obj%.*};        \
	done

clean:
	@- $(RM) $(PROGRAM_NAME) $(OBJECTS) *~ 

distclean: clean
