run: nogmp
	./nogmp

nogmp: nogmp.c arithmetic.c polynomial.c
	gcc -lgmp -o nogmp nogmp.c arithmetic.c polynomial.c

clean:
	rm nogmp

