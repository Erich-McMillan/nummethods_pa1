#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>


/* Function Declarations */
int rootBisection(double*, int, double, double, double*);
int rootNewton(double*, int, double, double, double*);
int Horners(double, double*, int, double*, double*);

int main(int _argc, char **_argv) {
	/*
		Description:
			Utilizes both Bisection and Newton-Raphson methods to find the root of a polynomial 	function within the specified bounds.

		Inputs:
			User will input the coefficients of the polynomial in monotonically declining order followed by the bounds in which the root is to be found, i.e. 2x^3 - x^1 + 3 = 0, with the bounds 0 and 1 shall be passed as ./program 2 0 -1 3 0 1
			_argv[1]: highest order coefficient
			_argv[.]:
			_argv[.]
			_argv[_argc-3]: lowest order coefficient
			_argv[_argc-2]: lower bound for search
			_argv[_argc-1]: upper bound for search

		Outputs:
			The results of the bounded Bisection search and Newton-Raphson search will be sent to stdout, the number of iterations for each search and the final result for the root as well as f(root) will be displayed next to the name of the search type. If an error occurs with either method the reason will be displayed. If either method fails to converge then the reason for such divergence will be stated.
	*/

	/* MARK: Variable Declarations */
	int polynomialorder=_argc-3, i, iterations_bisection, iterations_newton;
	double a, b, root_bisection, root_newton;
	double *polynomialValues;

	/* MARK: Allocate Memory */
	polynomialValues = malloc(polynomialorder * sizeof(double));

	/* MARK: Parse Commandline inputs */
	printf("Polynomial function: f(x) = ");
	for(i = 1; i < _argc-2; i++) {
		polynomialValues[i-1] = atof(_argv[i]);														// Get polynomial coefficients
		printf("\t%fx^%d +", polynomialValues[i-1], _argc-2-i-1);					// Print polynomial function
	}
	printf("\n");
	a = atof(_argv[_argc-2]);																						// Get a value
	b = atof(_argv[_argc-1]);																						// Get b value

	/* MARK: Bisection Method */
	iterations_bisection = rootBisection(polynomialValues, polynomialorder, a, b, &root_bisection);
	printf("\tBisection root of f(x) is %.9f. Which took %d iterations.\n", root_bisection, iterations_bisection);

	/* MARK: Newton's Method */
	iterations_newton = rootNewton(polynomialValues, polynomialorder, a, b, &root_newton);
	printf("\tNewton root of f(x) is %.9f. Which took %d iterations.\n", root_newton, iterations_newton);

	/* MARK: Clean-up and return */
	free(polynomialValues);
	return 0;
}

int rootBisection(double *_polyVals, int _polyorder, double _a, double _b, double *_root) {
	/*
		Description:
			Uses bisection to find the simple root of a polynomial function within the bounds of a and b. When f(a)*f(b) < 0 a root exists in the bounds [a,b]. The method can be broken into 4 steps:

			1. compute the midpoint between a and b: c. c = (b+a)/2 if |b-a|/2 is small enough (in our case we'd like to use the epsilon machine precision).
			2. If f(c)*f(a) < 0 then root is bounded between [a,c] let b->c
			3. Otherwise root must exist in [c,b] let a->c.
			4. Repeat continue to step 1. Program will end when epsilon is reached.

		Inputs:
			_polyVals (*double): the coefficients for the polynomial in monotonically decreasing order. The indicies are related to the order of coefficient where index 0 is the coefficient for the highest order and index n is the order of the nth polynomial order. If _polyVals is a NULL pointer then the function returns error code -1.

			_polyorder (int): the order of the polynomial function where _polyVals(_polyorder) returns the lowest order polynomial of x^0, and _polyVals(0) returns the highest order polynomial x^_polyorder-1. Therefore the order of the polynomial is _polyorder-1. If _polyorder is non-positive then the function returns an error code -2.


			_a (int): the lesser bound of the bisection search where _a < _b or the function returns an error code -3.

			_b (int): the greater bound of the bisection search where _b > _a or the function returns an error code -3.

			_root (*double): after the function is evaluated then the value for the root is placed at this memory address. If the function fails to find a root or otherwise fails the value at this address is unchanged. Function returns -4 if _root is a NULL pointer.

		Outputs:
			Upon success the function will return the number of iterations (always a positive number) needed to find the simple root of the polynomial and place the root at the memory address specified by _root. The error code below will be returned and an explaination for failure will be printed to STD_OUT.

			Error Codes*:
				-1: _polyVals is a NULL pointer
				-2: _polyorder is a non-positive integer
				-3: _a is not strictly less than or equal to _b
				-4: _root is a NULL pointer

				* note that more than one error code may return at a time. If more than one error occurs only the first will be caught and returned. The checks occur in the order presented above. If both _polyVals and _root are NULL pointers then only -1 will return.
	*/

	/* Check Inputs */
	if(_polyVals == NULL) {
		printf("\t~~in 'rootBisection' input _polyVals was NULL\n");
		return -1;
	}
	if(_polyorder <= 0) {
		printf("\t~~in 'rootBisection' input _polyorder was non-positive\n");
		return -2;
	}
	if(_a > _b) {
		printf("\t~~in 'rootBisection' input _a was not less than or equal to _b NULL\n");
		return -3;
	}
	if(_root == NULL) {
		printf("\t~~in 'rootBisection' input _root was NULL\n");
		return -4;
	}

	/* Declare variables */
	int iterations=0, ba;
	double curr_a = _a, curr_b = _b, curr_c, result_fa, result_fc;

	// try to compare (a + b)/2 == a or b

	/* Bisection to find root */
	//printf("a+b/2 = %f\ta = %f\tb=%f\n", (curr_a+curr_b)/2.0, curr_a, curr_b);
	while((curr_a+curr_b)/2.0 != curr_a && (curr_a+curr_b)/2.0 != curr_b) {	// Step 1. Check for epsilon
		curr_c = (curr_b+curr_a)/2.0;																			// Step 1. Compute midpoint
		Horners(curr_a ,_polyVals, _polyorder, &result_fa, NULL);	 				// Step 2. Calculate f(a)
		Horners(curr_c, _polyVals, _polyorder, &result_fc, NULL);					// Step 2. Calculate f(c)
		ba = (result_fa*result_fc) <= 0;
		if( ba ) {																												// Step 2. f(a)*f(c) < 0
			curr_b = curr_c;																								// Step 2. b->c
		} else {																													// Step 3. f(b)*f(c) < 0
			curr_a = curr_c;																								// Step 3. a->c
		}
		iterations++;																											// Step 4. Inc iter return Step 1.
		//printf("\ta = %f, b = %f, c = %f\n", curr_a, curr_b, curr_c);
	}

	/* Clean-up and return */
	*_root = curr_a;																										// Set _root to curr_a
	return iterations;
}
int rootNewton(double *_polyVals, int _polyorder, double _a, double _b, double *_root) {
	/*
		Description:
			Utilizes Newton-Raphson method for finding roots. This method requires that f(x) be differentiable and that you provide an initial guess _xO. The next guess is where f'(_xO) intersects the x-axis. The process is repeated until either the root is found or the search leaves the bounded area. The process for Newton-Raphson method is as follows:

				1. Start at some initial guess in our case we start at (a+b)/2
				2. Find the f(x) and f'(x) at this guess
				3. Next guess is nxt_x = curr_x - f(curr_x)/f'(curr_x)
				4. Repeat steps 2-4 until f(x) is within the machine precision

		Inputs:
			_polyVals (*double): the coefficients for the polynomial in monotonically decreasing order. The indicies are related to the order of coefficient where index 0 is the coefficient for the highest order and index n is the order of the nth polynomial order. If _polyVals is a NULL pointer then the function returns error code -1.

			_polyorder (int): the order of the polynomial function where _polyVals(_polyorder) returns the lowest order polynomial of x^0, and _polyVals(0) returns the highest order polynomial x^_polyorder-1. Therefore the order of the polynomial is _polyorder-1. If _polyorder is non-positive then the function returns an error code -2.

			_a (int): the lesser bound of the bisection search where _a < _b or the function returns an error code -3.

			_b (int): the greater bound of the bisection search where _b > _a or the function returns an error code -3.

			_root (*double): after the function is evaluated then the value for the root is placed at this memory address. If the function fails to find a root or otherwise fails the value at this address is unchanged. Function returns -4 if _root is a NULL pointer.

		Outputs:
			Upon success the function will return the number of iterations (always a positive number) needed to find the simple root of the polynomial and place the root at the memory address specified by _root. The error code below will be returned and an explaination for failure will be printed to STD_OUT.

			Error Codes*:
				-1: _polyVals is a NULL pointer
				-2: _polyorder is a non-positive integer
				-3: _a is not strictly less than or equal to _b
				-4: _root is a NULL pointer
				-5: if deriviative at some point is zero

				* note that more than one error code may return at a time. If more than one error occurs only the first will be caught and returned. The checks occur in the order presented above. If both _polyVals and _root are NULL pointers then only -1 will return.
	*/

	/* Check Inputs */
	if(_polyVals == NULL) {
		printf("\t~~in 'rootBisection' input _polyVals was NULL\n");
		return -1;
	}
	if(_polyorder <= 0) {
		printf("\t~~in 'rootBisection' input _polyorder was non-positive\n");
		return -2;
	}
	if(_a > _b) {
		printf("\t~~in 'rootBisection' input _a was not less than or equal to _b NULL\n");
		return -3;
	}
	if(_root == NULL) {
		printf("\t~~in 'rootBisection' input _root was NULL\n");
		return -4;
	}

	/* Declare variables */
	int numiter=0, bisection=0, biseciter=0, bounderr=0;
	double x_curr=(_a+_b)/2.0, x_prev, x_prevprev, fx, fdx, rt;

	/* Newton-Raphson algorithm */
	do {
		if(numiter > 1) {
			x_prevprev = x_prev;
		}
		x_prev = x_curr;
		Horners(x_curr, _polyVals, _polyorder, &fx, &fdx);								// Step 2. find f(x) f'(x)
		x_curr = x_curr - fx/fdx;																					// Step 3. find next guess
		numiter++;

		if(fdx == 0) {																										// Check fdx to see if min or max
			bounderr = 3;
			bisection = 1;
			if(numiter < 2) {
				biseciter = rootBisection(_polyVals, _polyorder, _a, _b, &rt);
			} else {
				biseciter = rootBisection(_polyVals, _polyorder, x_prev, x_prevprev, &rt);
			}
			continue;
		}
		if(x_curr == x_prevprev && numiter > 1) {													// Check for loop condition
			bisection = 1;
			if(x_curr < x_prev) {																						// Switch to bisection _xcurr <
				biseciter = rootBisection(_polyVals, _polyorder, x_curr, x_prev, &rt);
			} else {																												// Switch to bisection _xcurr >
				biseciter = rootBisection(_polyVals, _polyorder, x_prev, x_curr, &rt);
			}
		}
		if(x_curr > _b) {
			bounderr = 2;
		}
		if(x_curr < _a) {
			bounderr = 1;
		}
		//printf("\t\t%.18f\t%.18f\n", x_prev, x_curr);
	} while(x_curr != x_prev && !bisection && !bounderr);								// Step 4. check if converged

	/* Clean-up and return */
	switch(bounderr) {
		case 1:
			printf("\tNewton guess exceeded lower bound try new guess\n");
			break;
		case 2:
			printf("\tNewton guess exceeded upper bound try new guess\n");
			break;
		case 3:
			printf("\tNewton fdx was 0 attempted bisection\n");
			break;

	}
	if(bisection==0) {
		*_root = x_curr;
	} else {
		if(biseciter >= 0) {
			printf("\t\tNewton's utilized bisection with %d iterations\n", biseciter);
		}
		*_root = rt;
	}
	return numiter;
}
int Horners(double _x, double *_polyVals, int _polyorder, double *_fx, double *_fdx) {
	/*
		Description:
			Will evaluate the result of a polynomial function at specified value _x by carrying out synthetic division and determining the deriviative of the polynomial function at _x. The steps for determing the f(x) and f'(x) are:

			Step 1: Set p = _polyVals[_polyorder] and q = 0
   		Step 2: Do steps 3 and 4 for i from 0 to _polyorder-1, increasing by 1
      Step 3: set q = p + _x * q
      Step 4: set p = _polyVals[i] + _x * p
   		Step 5: The value of f(_x) is p and the value of f'(_x) is q

		Inputs:
			_x (double): the value at which the function will be evaluated

			_polyVals (*double): the coefficients for the polynomial in monotonically decreasing order. The indicies are related to the order of coefficient where index 0 is the coefficient for the highest order and index n is the order of the nth polynomial order. If _polyVals is a NULL pointer then the function returns error code -1.

			_polyorder (int): the order of the polynomial function where _polyVals(_polyorder) returns the lowest order polynomial of x^0, and _polyVals(0) returns the highest order polynomial x^_polyorder-1. If _polyorder is non-positive then the function returns an error code -2.

			_fx (*double): A pointer to where the result of f(x) will be stored. If the function fails due to a programmed error code _fx is not changed. If _fx and _fdx are NULL error code -3.

			_fdx (*double): A pointer to where the result of f'(x) will be stored. If the function fails due to a programmed error code _fdx is not changed. If _fx and _fdx are NULL error code -3.

		Outputs:
			Places the f(_x) -> *_fx and f'(_x) -> *_fdx. By passing a NULL value for a single input that input will not be returned. Returns 0 if all is sucessful and the following error codes if any issues occur.

			Error Codes:
				-1: _polyVals is a NULL pointer
				-2: _polyorder is a non-positive integer
				-3:	_fx is a NULL pointer AND _fdx is a NULL pointer
	*/

	/* Check Inputs */
	if(_polyVals == NULL) {
		printf("\t~~in 'rootBisection' input _polyVals was NULL\n");
		return -1;
	}
	if(_polyorder <= 0) {
		printf("\t~~in 'rootBisection' input _polyorder was non-positive\n");
		return -2;
	}
	if(_fx == NULL && _fdx == NULL) {
		printf("\t~~in 'rootBisection' input _fx and _fdx are BOTH NULL\n");
		return -3;
	}

	/* Declare Variables */
	int i = 0;
	double p=_polyVals[0], q=0;																					// Step 1. q=0, p=_polyVals[order]


	/* Horner's method */
	for(i = 1; i < _polyorder; i++) {																		// Step 2. i=0 -> _polyorder
		q = p + _x * q;																										// Step 3: q=p+_x*q
		p = _polyVals[i] + _x*p;																					// Step 4: p=_polyVals[i]+_x*p
	}

	/* Clean-up and return */
	if(_fx != NULL) *_fx = p;
	if(_fdx != NULL) *_fdx = q;
	return 0;
}
