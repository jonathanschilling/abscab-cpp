
CC?=gcc
CFLAGS=--std=c99

ABSCAB_DIR=abscab
TEST_DIR=test

.PHONY: all test clean demo_programs test_programs

all: test_programs demo_programs

test_programs: test_cel test_abscab

demo_programs: demo_abscab

test: test_programs
	@ ./test_cel
	@ ./test_abscab

demo: demo_programs
	@ ./demo_abscab

clean:
	rm -f test_cel test_abscab demo_abscab

test_cel: $(ABSCAB_DIR)/cel.h $(TEST_DIR)/util.h $(TEST_DIR)/test_cel.c
	$(CC) $(CFLAGS) $(TEST_DIR)/test_cel.c -I$(ABSCAB_DIR) -o test_cel -lm

test_abscab: $(ABSCAB_DIR)/abscab.h $(TEST_DIR)/util.h $(TEST_DIR)/test_abscab.c
	$(CC) $(CFLAGS) $(TEST_DIR)/test_abscab.c -I$(ABSCAB_DIR) -o test_abscab -lm

demo_abscab: $(ABSCAB_DIR)/abscab.h $(TEST_DIR)/util.h $(TEST_DIR)/demo_abscab.c
	$(CC) $(CFLAGS) $(TEST_DIR)/demo_abscab.c -I$(ABSCAB_DIR) -o demo_abscab -lm -O3 -fopenmp
