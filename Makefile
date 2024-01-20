run: aks
	./aks

aks: aks.c arithmetic.c polynomial.c
	gcc -lm -lgmp -o aks aks.c arithmetic.c polynomial.c

clean:
	rm aks

