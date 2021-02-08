#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>

typedef struct Ratnums {
	   short upper;
	   short lower;
	} ratnum;


//Performs Fourier–Motzkin elimination, entry into a programming competitions scored by execution speed
bool fm(size_t rows, size_t cols, signed char a[rows][cols], signed char c[rows])
{
	bool realloc = false;
	void free_all(ratnum *q, ratnum *qprim, ratnum **t, ratnum **tprim, int maxs){
		free(t);
		free(q);
		free(tprim);
		free(qprim);
	}
	bool altint(ratnum a, int b){
		bool ret = a.upper < b*a.lower;
		if(a.lower < 0 && a.upper != 0){
			ret = !ret;
			
		}
		return ret;
	}
	bool aeqint(ratnum a, int b){
		bool ret = a.upper == b*a.lower;
		return ret;
	}
	bool agtint(ratnum a, int b){
		bool ret = a.upper > b*a.lower;
		if(a.lower < 0 && a.upper != 0){
			ret = !ret;
			
		}
		return ret;
	}

	bool altb(ratnum a, ratnum b){
		bool ret = a.upper*b.lower < b.upper*a.lower;
		if((a.lower < 0 || b.lower < 0) && !(a.lower < 0 && b.lower < 0)){
			ret = !ret;
		}
		return ret;
	}


	ratnum subtract(ratnum one, ratnum two){
		short a = one.upper;
		short b = two.upper;
		short c = one.lower;
		short d = two.lower;
		long up = a*d - b*c;
		long down = c*d;
		while(up > 32767 || up < -32767 || down > 32767 || down < -32767){
			up = up >> 1;
			down = down >> 1;
		}
		one.upper = (short)up;
		one.lower = (short)down;
		return one;
	}

	ratnum divide(ratnum one, ratnum two){
		short a = one.upper;
		short b = two.upper;
		short c = one.lower;
		short d = two.lower;
		long up = a*d;
		long down = b*c;
		while(up > 32767 || up < -32767 || down > 32767 || down < -32767){
			up = up >> 1;
			down = down >> 1;
		}
		one.upper = (short)up;
		one.lower = (short)down;
		return one;
	}




	//This is a guess for the largest value of s
	int maxs = rows*1.4;

	//Inializes the matrices
	ratnum *q = (ratnum*)malloc(maxs*sizeof(ratnum));
	ratnum *qprim = (ratnum*)malloc(maxs*sizeof(ratnum));

	ratnum *t = (ratnum *)malloc(maxs * cols * sizeof(ratnum));
	ratnum *tprim = (ratnum *)malloc(maxs * cols * sizeof(ratnum));


	//Transfers to input to swap-matrices
	for(int j = 0; j < rows; j++){		
		for(int i = 0; i < cols; i++){
			tprim[j*cols+i].upper = (short)a[j][i];
			tprim[j*cols+i].lower = 1;
		}
		qprim[j].upper = (short)c[j];
		qprim[j].lower = 1;
	}
	
	//Initial value for s
	int s = rows;



	//Main loop
	for (int r = cols; r > 0; r--){

		//Transfers from swap-matrix to main-matrix i order and counts (+,-,0)
		int n = 0;
		for(int i = 0; i < s; i++){
			if(agtint(tprim[i*cols+r-1], 0)){
				for(int j = 0; j < r; j++){
					t[n*cols+j] = tprim[i*cols+j];
					q[n] = qprim[i];
				}
				n++;
			}
		}
		int n1 = n;
		for(int i = 0; i < s; i++){
			if(altint(tprim[i*cols+r-1], 0)){
				for(int j = 0; j < r; j++){
					t[n*cols+j] = tprim[i*cols+j];
					q[n] = qprim[i];
				}
				n++;
			}
		}
		int n2 = n;

		for(int i = 0; i < s; i++){
			if(aeqint(tprim[i*cols+r-1], 0)){
				for(int j = 0; j < r; j++){
					t[n*cols+j] = tprim[i*cols+j];
					q[n] = qprim[i];
				}
				n++;
			}
		}	
		//Normalizes x_r
		for(int i = 0; i < s; i++){
			for(int j = 0; j < r-1; j++){	
				t[i*cols+j] = divide(t[i*cols+j], t[i*cols+r-1]);
			}
			q[i] = divide(q[i], t[i*cols+r-1]);
		}

		//We done
		if(r == 1){
			for(int i = 0; i < n1; i++){		
				for(int j = n1; j < n2; j++){
					if(altb(q[i], q[j])){
						free_all(q, qprim, t, tprim, maxs);
						return 0;
					}
				}
			}
			for(int i = n2; i < s; i++){
				if(altint(q[i], 0)){
					free_all(q, qprim, t, tprim, maxs);
					return 0;
				}
			}
			free_all(q, qprim, t, tprim, maxs);
			return 1;
		}

		//Number of equations in next round
		int sprim = s - n2 + n1*(n2-n1);
		if(sprim > maxs){
			realloc = true;
			maxs = sprim;
			free(qprim);
			free(tprim);
			qprim = (ratnum*)malloc(maxs*sizeof(ratnum));
			tprim = (ratnum*)malloc(maxs * cols * sizeof(ratnum));
		}
		//No equations = no problems
		if(sprim == 0) {
			free_all(q, qprim, t, tprim, maxs);
			return 1;
		}

		//Create the new equations
		int m = 0;
		for(int i = 0; i<n1; i++){
			for(int j = n1; j < n2; j++){
				qprim[m] = subtract(q[i],q[j]);
				for(int z = 0; z < r; z++){
					tprim[m*cols+z] = subtract(t[j*cols+z], t[i*cols+z]); 
				}
				m++;
			}
		}
		for(int i = n2; i<s; i++){
			for(int z = 0; z < r; z++){
				tprim[m*cols+z] = t[i*cols+z]; 
			}
			m++;
		}
		s=sprim;
		//If we underestimated the size of s we need to allocate more memory
		if(realloc){
			free(q);
			free(t);
			q = (ratnum*)malloc(maxs*sizeof(ratnum));
			t = (ratnum *)malloc(maxs * cols * sizeof(ratnum));
			realloc = false;
		}
	}
}
