ifeq ($(CC),)
CC = gcc
endif

CFLAGS = -Wall -g

THIS_DIR = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

BIN = bin
OBJ = obj
SRC = src
TEST = test
TEST_BIN = $(TEST)/bin
VENDOR = vendor
UNITY = $(VENDOR)/unity/src
TOMMY = $(VENDOR)/tommyds

OBJS := $(OBJ)/array.o

.PHONY: clean apple_valgrind

apple:
	$(CC) $(CFLAGS) -o $(BIN)/$@ $(SRC)/$@.c

apple_docker:
	docker run --workdir $(HOME) --entrypoint $(CC) -v $(THIS_DIR):$(HOME) mooreryan/valgrind $(CFLAGS) -o $(BIN)/$@ $(SRC)/apple.c

apple_valgrind: apple_docker
	docker run --workdir $(HOME) -v $(THIS_DIR):$(HOME) mooreryan/valgrind $(BIN)/$^

$(OBJ)/array.o:
	gcc -o $(OBJ)/array.o -c $(SRC)/array.c

$(OBJ)/file.o:
	gcc -Ivendor/tommyds/ -o $(OBJ)/file.o -c $(SRC)/file.c

$(UNITY)/unity.o:
	gcc -o $(UNITY)/unity.o -c $(UNITY)/unity.c

$(OBJ)/tommyarray.o:
	gcc -o $(OBJ)/tommyarray.o -c $(TOMMY)/tommyarray.c

test_array: $(OBJ)/array.o $(UNITY)/unity.o
	gcc -Ivendor/ruby_like_c/ -Ivendor/tommyds/ -Ivendor/unity/src -Isrc -o FFFF test/test_array.c $^
	./FFFF

test_file: $(OBJ)/file.o $(UNITY)/unity.o $(OBJ)/tommyarray.o
	gcc -Ivendor/ruby_like_c/ -Ivendor/tommyds/ -Ivendor/unity/src -Isrc -o FFFF test/test_file.c $^
	./FFFF

test_file_valgrind:
	docker run --workdir $(HOME) --entrypoint $(CC) -v $(THIS_DIR):$(HOME) mooreryan/valgrind $(CFLAGS) -Ivendor/tommyds/ -o $(OBJ)/file.o -c $(SRC)/file.c
	docker run --workdir $(HOME) --entrypoint $(CC) -v $(THIS_DIR):$(HOME) mooreryan/valgrind $(CFLAGS) -o $(UNITY)/unity.o -c $(UNITY)/unity.c
	docker run --workdir $(HOME) --entrypoint $(CC) -v $(THIS_DIR):$(HOME) mooreryan/valgrind $(CFLAGS) -o $(OBJ)/tommyarray.o -c $(TOMMY)/tommyarray.c
	docker run --workdir $(HOME) --entrypoint $(CC) -v $(THIS_DIR):$(HOME) mooreryan/valgrind $(CFLAGS) -Ivendor/ruby_like_c/ -Ivendor/tommyds/ -Ivendor/unity/src -Isrc -o FFFF test/test_file.c $(OBJ)/file.o $(UNITY)/unity.o $(OBJ)/tommyarray.o
	docker run --workdir $(HOME) -v $(THIS_DIR):$(HOME) mooreryan/valgrind --leak-check=full ./FFFF

clean:
	rm -r $(BIN)/* $(OBJ)/* $(UNITY)/unity.o
