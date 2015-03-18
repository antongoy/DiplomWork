all: ranlib.o main.o com.o
	gcc -Wall com.o ranlib.o main.o -o residual_min -llapacke -lm -llapack

ranlib.o: ranlib/ranlib.c
	gcc -O3 -c ranlib/ranlib.c 

com.o: ranlib/com.c
	gcc -O3 -c ranlib/com.c 

main.o: main.c
	gcc -Wall -g -ggdb -c main.c

clean:
	rm -f *.o


