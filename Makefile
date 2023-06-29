# testing
TARGET = clara
SRC = ./src
INC = ./include

CC_STANDARD = -std=c++11
CC = g++
WARNINGS = -Wall -Wextra
MULTIPROC = -fopenmp
OPTIM = -mtune=native -msse2
EIGEN = /usr/include/eigen3

CFLAGS = -c -pedantic $(CC_STANDARD) $(WARNINGS) $(MULTIPROC) $(OPTIM)\
				 -isystem $(EIGEN) -I $(INC)

CFLAGS_RELEASE = -O2 -DNDEBUG -DEIGEN_NO_DEBUG # Release flags
CFLAGS_DEBUG = -DDEBUG -g3 # Debug flags

# Use gomp multi-processing library
LDFLAGS = -lgomp

SOURCES = $(wildcard $(SRC)/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: release

# Release make configuration
release: CFLAGS += $(CFLAGS_RELEASE)
release: LDFLAGS += $(MULTIPROC)
release: $(SOURCES) $(TARGET)

# Debug make configuration
debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LDFLAGS += $(MULTIPROC)
debug: $(SOURCES) $(TARGET) 


$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

# Clean-up
clean:
	@echo 'Removing...'
	@rm -fv $(SRC)/*.o 
	@rm -fv $(TARGET)
