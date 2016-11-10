SRCS = \
			 main.cpp \
			 blas.cpp \
			 io.cpp \
			 selecter.cpp \
			 start.cpp \
			 tools.cpp \
			 CRS/cg.cpp \
			 CRS/cr.cpp \
			 CRS/gcr.cpp \
			 CRS/gmres.cpp \
			 CRS/kskipcg.cpp \
			 CRS/vpcg.cpp \
       CRS/vpcr.cpp \
       CRS/vpgcr.cpp \
       CRS/vpgmres.cpp \
       CRS/bicg.cpp \
       CRS/kskipbicg.cpp

# CRS/kskipcr.c \

INC_DIR = .

BUILD_DIR = Build

OBJS=$(addprefix $(BUILD_DIR)/,$(patsubst %.cpp,%.o,$(SRCS)))
DEPS=$(patsubst %.o,%.d, $(OBJS))

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	# CC = gcc
	CC = g++
	# CFLAGS += -lm -std=c99 -pg
	CFLAGS += -lm
endif
ifeq ($(UNAME), Darwin)
	CC = gcc-6
endif
ifeq ($(UNAME), Windows_NT)
	exit 0
endif

TARGET = Solver
CFLAGS += $(addprefix -I,$(INC_DIR)) -fopenmp -Wall -fno-use-linker-plugin
LDFLAGS +=

.PHONY: all clean cpu debug

all: cpu

cpu: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

$(BUILD_DIR)/%.o : %.cpp
	mkdir -p $(dir $@); \
		$(CC) -c $(CFLAGS) -o $@ $<

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

debug:
	gprof ./$(TARGET) | ./gprof2dot.py | dot -Tpng -o debug.png

.PHONY: all clean
-include $(DEPS)
