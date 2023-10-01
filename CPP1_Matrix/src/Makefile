

STD_FLAGS=-Wall -Wextra -Werror -std=c++17 -O0 -g 
DEBUG_FLAGS=-O0 -g 
TEST_FLAGS=-lgtest -lgtest_main -pthread
CC=gcc
LIB=s21_matrix_oop.a
CPP_LIB=-lstdc++

all: clean test

$(LIB): s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o

test: clean $(LIB) test.o
	$(CC) $(STD_FLAGS) -o test.out  test.o $(LIB) $(TEST_FLAGS) $(CPP_LIB)
	./test.out

clean:
	rm -rf *.a *.out *.o


# Вспомогательные цели

memcheck:  clean $(LIB) test.o 
	$(CC) $(STD_FLAGS) -o test.out  test.o $(LIB) $(TEST_FLAGS) $(CPP_LIB)
	valgrind --tool=memcheck --leak-check=yes  ./test.out

stylecheck:
	clang-format -i *.cpp *.h
	clang-format -n *.cpp *.h 

s21_matrix_oop.o:
	$(CC)  $(STD_FLAGS) -c s21_matrix_oop.cpp $(CPP_LIB)

test.o: 
	$(CC)  $(STD_FLAGS) -c test.cpp $(CPP_LIB)
