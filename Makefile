run: aks
	./aks

aks: aks.c arithmetic.c polynomial.c
	gcc -lgmp -o aks aks.c arithmetic.c polynomial.c

clean:
	rm nogmp

