all: ranlib.o main.o com.o
	gcc -Wall -O3 com.o ranlib.o main.o -o residual_min -llapacke -lm -llapack

ranlib.o: ranlib/ranlib.c
	gcc -O3 -c ranlib/ranlib.c 

com.o: ranlib/com.c
	gcc -O3 -c ranlib/com.c

$FLAGS = -g -ggdb -O0

main.o: main.c
	gcc -Wall -O3  -c main.c

clean:
	rm -f *.o


