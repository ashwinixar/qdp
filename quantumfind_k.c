/*
*****************************************************
*	About: The code implements Quantum algorithm 	*
*	to find division property of a multiset			*
*	Usage: Run in command prompt "quantumfind_k.exe"*
*	Author: Ashwini Kumar Malviya					*
*	Email: ashwinixar@gmail.com						*
*****************************************************
*/


#include <stdio.h>
#include <math.h>
#include "QuEST.h"

//Multiset
const int X_size = 6;
int multiset_X[6] = { 0x1, 0x1, 0x1, 0x3, 0x3, 0x5 };

//Evaluates bit product of "x" and "u" i.e. pi_u(x)
int bit_product(int x, int u)
{
	int x_u = 1;
	while(u)
	{
		if(u & 1) x_u &= (x & 1);
		u >>= 1;
		x >>= 1;
	}
	return x_u;
}

//Find hamming weight of a given number
int ham_wgt(int x)
{
	int wgt = 0;
	while(x)
	{
		if(x & 1) wgt++;
		x >>= 1;
	}
	return wgt;
}

//Finds the next number (less than n) after "x" with a given hamming weight "w"
//If no number is found then 0 is returned
int next_num(int n, int x, int w)
{
	for(int i = x + 1; i < n; i++)
		if(ham_wgt(i) == w) return i;
	return 0;
}

double count(int u)
{
	QuESTEnv env = createQuESTEnv();

    //Value of Pi
    const double pi = 3.14;

    const int p = 6, P = 64;
    const int n = 3, N = 8;
    Qureg qubits = createQureg((p + n), env);
    initZeroState(qubits);
    //reportQuregParams(qubits);
    //reportQuESTEnv(env);

    //Unitary matrix creation
	ComplexMatrixN e = createComplexMatrixN(n);
    for(int i = 0; i < 8; i++)
    {
    	if(i >= X_size) e.real[i][i] = 1;
    	else
    	{
    		if(bit_product(multiset_X[i], u)) e.real[i][i] = -1;
    	    else e.real[i][i] = 1;
    	}
    }

    //Apply Hadamard on all the (p + n)-qubits to create superposition
    //Here each element in last n-qubit superposition is considered to represent "x" 
    for(int i = 0; i < (p + n); i++)
        hadamard(qubits, i);
    
    //Target qubits on which grover's search must be performed to count the number of correct result
    int targs[n];
    for(int i = 0; i < n; i++)
    	targs[i] = p + i;

    //Applying controlled grover's operator
    int iterations = 1;
    for(int gi = p - 1; gi >= 0; gi--)
    {
    	int ctrls[] = { gi };
    	for(int j = 0; j < iterations; j++)
    	{
    		//Marking
	    	multiControlledMultiQubitUnitary(qubits, ctrls, 1, targs, n, e);
	    	
	        //Diffusion
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        multiControlledPhaseFlip(qubits, targs, n);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
    	}
    	iterations *= 2;
    }
	
    //Inverse QFT on the first p-qubits
    for(int i = 0; i < p; i++)
    {
    	for(int j = 0; j < i; j++)
    	{
    		qreal angle = -2.0 * pi / (pow(2, i - j + 1));
    		controlledPhaseShift(qubits, j, i, angle);
    	}
    	hadamard(qubits, i);
    }
    
    //Measures the qubits to get the value of first p-qubits
    qreal prob = 0, max_prob = 0;
    int p_val = 0;
    for(int i = 0; i < P; i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob > max_prob)
    	{
    		max_prob = prob;
    		p_val = i;
    	}
    }
    if(p_val > (P / 2)) p_val = P - p_val;
    //printf("\nThe measured value of (f_tilde) is %d", p_val);
    //printf("\nThe approximate count of the correct answers is: %f\n", (N * pow(sin(p_val * pi / P), 2)));

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return (N * pow(sin(p_val * pi / P), 2));
}

int main (int narg, char *varg[])
{
    int k;
    int found = 0;
    int u = 0;
    int num = 8; //2^3
    for(int i = 1; i <= 3; i++) //Running for all hamming weights
    {
    	k = i;
    	u = next_num(num, u, i);
    	while(u)
    	{
		    if((int)round(count(u)) % 2)
		    {
		    	found = 1;
		    	break;
		    }
		    if(found) break;
			u = next_num(num, u, i);
    	}
    	if(found) break;
    	u = 0;
    }
    
    if(k == 1 && X_size % 2 == 1) k = 0;
    printf("\nThe Division Property holds for k = %d when u = %02X\n", k, u);

    return 0;
}
