.PHONY: all lsolver rebuild clean testall run_test% dist

GPP ?= g++
OTHERFLAGS+=
IDIR=include
SELDON=lib/seldon
override CFLAGS+=-I$(IDIR) -I$(SELDON) -g -Wall $(OTHERFLAGS) -std=c++1z -fopenmp -Wno-format-extra-args
PEDANTIC_CFLAGS=-std=c++1z -Werror -Wpedantic -Wall -Wextra -Wformat=2 -O -Wuninitialized -Winit-self -Wswitch-enum -Wdeclaration-after-statement -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Wno-long-long

ODIR=obj
LDIR =lib
SRCDIR=src
TSTDIR=test

LDFLAGS=-lm

LIBS = $(wildcard $(LDIR)/*.o)
DEPS = $(wildcard $(IDIR)/*.h)
SRCS = $(wildcard $(SRCDIR)/*.cpp)
_OBJ = $(patsubst $(SRCDIR)/%,$(ODIR)/%,$(SRCS:.cpp=.o))
OBJ = $(filter-out $(ODIR)/run.o,$(_OBJ))

DEPS2 := $(OBJ:.o=.d)
-include $(DEPS2)

all: lsolver

$(ODIR)/%.o: $(SRCDIR)/%.cpp
	+@[ -d $(ODIR) ] || mkdir -p $(ODIR)
	$(GPP) -MMD $(CFLAGS) -c -o $@ $< 

lsolver: $(OBJ) $(ODIR)/run.o
	echo $(SRCS)
	$(GPP) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LIBS) 

clean:
	-rm -f $(ODIR)/*.o *~ core.* $(INCDIR)/*~
	-rm -f $(ODIR)/*.d
	-rm -f lsolver
	-rm -f test1 test2 test3 test4 test5 test6
	-rm -f dist.tar.gz
	-rm -rf profdata/


dist: clean
	tar -cvzf dist.tar.gz src/* include/* Makefile README.md

test%: $(OBJ) $(TSTDIR)/test%.cpp
	$(GPP) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

run_test%: test%
	./$<

rebuild: clean all

testall: clean run_test1 run_test2 run_test3 run_test4 run_test5 run_test6