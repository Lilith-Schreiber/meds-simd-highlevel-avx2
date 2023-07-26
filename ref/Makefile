SHELL := /bin/bash

# PARAM = MEDS167717
PARAM = MEDS201254
# PARAM = MEDS134180
# PARAM = toy

LIBS = -lssl -lcrypto
CC := gcc

ifdef DEBUG
CFLAGS := -g -Wall -D$(PARAM) -DDEBUG -mavx512f
OBJDIR := debug/$(PARAM)
else
CFLAGS := -O3 -Wall -D$(PARAM) -mavx512f
OBJDIR := build/$(PARAM)
endif

ifdef LOG
CFLAGS += -DDEBUG
endif

CFLAGS += -I$(OBJDIR) -INIST

EXES = test bench KAT_test

TARGETS := ${EXES:%=$(OBJDIR)/%}

.PHONY: default all clean

default: $(EXES)

# OBJECTS = meds.o util.o seed.o osfreq.o fips202.o matrixmod.o bitstream.o randombytes.o
OBJECTS = meds.o util.o seed.o fips202.o matrixmod.o bitstream.o randombytes.o measure.o vec.o
HEADERS = $(wildcard *.h) params.h $(OBJDIR)/api.h

BUILDOBJ := ${OBJECTS:%=$(OBJDIR)/%}


params.h: params.py
	python $< -p > $@

$(OBJDIR)/api.h: params.py $(OBJDIR)
	python $< -a $(PARAM) > $@


$(EXES) : % :
	@make $(OBJDIR)/$(@F)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BUILDOBJ) : $(OBJDIR)/%.o: %.c $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGETS) : $(OBJDIR)/%: %.c $(BUILDOBJ)
	$(CC) $(@F).c $(BUILDOBJ) $(CFLAGS) $(LIBS) -o $@

RUN: test
	@echo ""
	@$(OBJDIR)/test

RUN_ALL:
	for p in `./params.py -l`; do \
	echo "Running par set" $$p; \
	make RUN PARAM=$$p; \
	echo ""; \
	done;

.PHONY: debug KAT

KAT: 
	$(MAKE) KAT_test PARAM=$(PARAM)
	build/$(PARAM)/KAT_test

KAT_BUILD: 
	$(MAKE) KAT_test PARAM=$(PARAM)

KAT_RUN: 
	build/$(PARAM)/KAT_test

KAT_BUILD_ALL:
	@for p in `./params.py -l`; do \
	echo "Building par set" $$p; \
	make PARAM=$$p KAT_BUILD; \
	echo ""; \
	done

KAT_RUN_ALL:
	@for p in `./params.py -l`; do \
	echo "Running par set" $$p; \
	make PARAM=$$p KAT_RUN; \
	echo ""; \
	done

KAT_ALL:
	make KAT_BUILD_ALL
	make KAT_RUN_ALL

BENCH:
	make bench PARAM=$(PARAM)
	$(OBJDIR)/bench | ./proc_bench.py >> bench.txt

BENCH_ALL:
	rm -f bench.txt
	for p in `./params.py -l`; do \
	echo "Running par set" $$p; \
	make PARAM=$$p BENCH; \
	echo ""; \
	cat bench.txt; \
	done

all:
	for p in `./params.py -l`; do \
	make PARAM=$$p; \
	make PARAM=$$p DEBUG=true; \
	done


clean:
	rm -rf build/ debug/ PQCsignKAT_*

