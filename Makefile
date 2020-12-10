# ------------------------------------------------------------------------------
# default settings

ifndef USE_OMP
USE_OMP = 1
endif

ifndef INDEX_LENGTH
INDEX_LENGTH = U32
endif

ifndef TIMING
TIMING = 0
endif

# set zero when the machine does not support SSE 4.0
WITH_STNI = 1

# ------------------------------------------------------------------------------
# compiler flags
CC	= gcc
CFLAGS	= -Wall -O3 -std=c99 -ftree-vectorize -msse4.2 -funroll-loops
LDFLAGS	= -Llib/
LIBS	= -lgomp -lm

ifeq ($(USE_OMP), 1)
CFLAGS	+= -fopenmp -D__OPENMP__
SUF_OMP	= _omp
endif

ifeq ($(TIMING), 1)
CFLAGS  += -D__TIMING__
endif

ifeq ($(WITH_STNI), 1)
CFLAGS  += -D__WITH_STNI__
endif

ifeq ($(INDEX_LENGTH), U32)
CFLAGS	+= -D__INDEX_U32__
SUF_KEY	= _u32
endif

ifeq ($(INDEX_LENGTH), U64)
CFLAGS	+= -D__INDEX_U64__
SUF_KEY	= _u64
endif


# ------------------------------------------------------------------------------
# definitions 
EXE	= pasqual$(SUF_OMP)$(SUF_KEY)
SRC     = src
OBJ     = obj
BIN	= bin
INS	= -I$(SRC)
FULLSOURCES	= $(wildcard $(SRC)/*.c)
SOURCES		= $(notdir $(FULLSOURCES))
OBJS		= $(patsubst %.c, %.o, $(SOURCES))
OBJECTS		= $(addprefix $(OBJ)/, $(OBJS))


# ------------------------------------------------------------------------------
# targets
all:			$(EXE)

$(EXE):			$(BIN)/$(EXE)

# rules for executables
$(BIN)/$(EXE): $(OBJECTS)
	$(CC) $(CFLAGS) $(LIBS) $(INS) $(LDFLAGS) -o $@ $^

# rules for source files
$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) $(INS) -o $@ -c $<

clean:
	rm -f $(BIN)/$(EXE) $(OBJECTS)
