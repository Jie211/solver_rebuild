SRCS = \
			 main.c \
			 blas.c \
			 io.c \
			 selecter.c \
			 start.c \
			 tools.c \
			 CRS/cg.c \
			 CRS/cr.c \
			 CRS/gcr.c \
			 CRS/gmres.c \
			 CRS/kskipcg.c \
			 CRS/kskipcr.c \
			 CRS/vpcg.c \
       CRS/vpcr.c \
       CRS/vpgcr.c \
       CRS/vpgmres.c \
       CRS/bicg.c

INC_DIR = .

BUILD_DIR = Build

OBJS=$(addprefix $(BUILD_DIR)/,$(patsubst %.c,%.o,$(SRCS)))
DEPS=$(patsubst %.o,%.d, $(OBJS))

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	CC = gcc
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

$(BUILD_DIR)/%.o : %.c
	mkdir -p $(dir $@); \
		$(CC) -c $(CFLAGS) -o $@ $<

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean
-include $(DEPS)
