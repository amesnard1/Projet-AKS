run: nogmp
	./nogmp

nogmp: nogmp.c arithmetic.c polynomial.c
	gcc -lm -o nogmp nogmp.c arithmetic.c polynomial.c

clean:
	rm nogmp

