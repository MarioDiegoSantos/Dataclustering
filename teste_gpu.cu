#include <cuda.h>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <sstream>
#include "gpu_cgrasp.h"


/********* XINSHEYANG2 ************/

__device__ float d_XINSHEYANG2(float *x, int n) {
	float sum1 = 0;
	for(int i = 0; i < n; i++) {
		sum1 += fabs(x[i]);
	}
	
	float sum2 = 0;
	for(int i = 0; i < n; i++) {
		sum2 += sinf(powf(x[i],2));
	}

	return sum1 * expf(-sum2);
}

float h_XINSHEYANG2(float *x, int n) {
	float sum1 = 0;
	for(int i = 0; i < n; i++) {
		sum1 += fabs(x[i]);
	}
	
	float sum2 = 0;
	for(int i = 0; i < n; i++) {
		sum2 += sinf(powf(x[i],2));
	}

	return sum1 * expf(-sum2);
}

__device__ fptr_t d_p_XINSHEYANG2 = d_XINSHEYANG2;

/********* PINTER ************/

__device__ float d_PINTER(float *x, int n) {
	float sum1 = 0, sum2 = 0, sum3 = 0;

	for(int i = 0; i < n; i++) {
		float A, B;
		if(i==0) {
			A = x[n-1] * sinf(x[i]) + sinf(x[i+1]);
			B = powf(x[n-1], 2) - 2*x[i] + 3*x[i+1] - cosf(x[i]) + 1;
		} else if(i==n-1) {
			A = x[i-1] * sinf(x[i]) + sinf(x[0]);
			B = powf(x[i-1], 2) - 2*x[i] + 3*x[0] - cosf(x[i]) + 1;
		} else {
			A = x[i-1] * sinf(x[i]) + sinf(x[i+1]);
			B = powf(x[i-1], 2) - 2*x[i] + 3*x[i+1] - cosf(x[i]) + 1;
		}
		
		sum1 += (i+1) * powf(x[i],2);
		sum2 += 20 * (i+1) * powf(sinf(A),2);
		sum3 += (i+1) * log10f(1+(i+1)*powf(B,2));
	}

	return sum1 + sum2 + sum3;
}

float h_PINTER(float *x, int n) {
	float sum1 = 0, sum2 = 0, sum3 = 0;

	for(int i = 0; i < n; i++) {
		float A, B;
		if(i==0) {
			A = x[n-1] * sinf(x[i]) + sinf(x[i+1]);
			B = powf(x[n-1], 2) - 2*x[i] + 3*x[i+1] - cosf(x[i]) + 1;
		} else if(i==n-1) {
			A = x[i-1] * sinf(x[i]) + sinf(x[0]);
			B = powf(x[i-1], 2) - 2*x[i] + 3*x[0] - cosf(x[i]) + 1;
		} else {
			A = x[i-1] * sinf(x[i]) + sinf(x[i+1]);
			B = powf(x[i-1], 2) - 2*x[i] + 3*x[i+1] - cosf(x[i]) + 1;
		}
		
		sum1 += (i+1) * powf(x[i],2);
		sum2 += 20 * (i+1) * powf(sinf(A),2);
		sum3 += (i+1) * log10f(1+(i+1)*powf(B,2));
	}

	return sum1 + sum2 + sum3;
}

__device__ fptr_t d_p_PINTER = d_PINTER;

/********* PATHOLOGICAL ************/

__device__ float d_PATHOLOGICAL(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		float num = powf(sinf(sqrtf(100* powf(x[i],2) + powf(x[i+1],2))),2) - 0.5;
		float den = 1 + 0.001*powf(powf(x[i],2)-2*x[i]*x[i+1]+powf(x[i+1],2),2);
		sum += 0.5 + num/den;
	}
	return sum;
}

float h_PATHOLOGICAL(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		float num = powf(sinf(sqrtf(100* powf(x[i],2) + powf(x[i+1],2))),2) - 0.5;
		float den = 1 + 0.001*powf(powf(x[i],2)-2*x[i]*x[i+1]+powf(x[i+1],2),2);
		sum += 0.5 + num/den;
	}
	return sum;
}

__device__ fptr_t d_p_PATHOLOGICAL = d_PATHOLOGICAL;

/********* ZAKHAROV ************/

__device__ float d_ZAKHAROV(float *x, int n) {
	float sum1 = 0;
	for(int i = 0; i < n; i++) {
		sum1 += powf(x[i],2);
	}

	float sum2 = 0;
	float sum3 = 0;
	for(int i = 0; i < n; i++) {
		sum2 += (i+1)*x[i];
	}
	sum2 = 1/2 * sum2;
	sum3 = powf(sum2,4);
	sum2 = powf(sum2,2);

	return sum1+sum2+sum3;
}

float h_ZAKHAROV(float *x, int n) {
	float sum1 = 0;
	for(int i = 0; i < n; i++) {
		sum1 += powf(x[i],2);
	}

	float sum2 = 0;
	float sum3 = 0;
	for(int i = 0; i < n; i++) {
		sum2 += (i+1)*x[i];
	}
	sum2 = 1/2 * sum2;
	sum3 = powf(sum2,4);
	sum2 = powf(sum2,2);

	return sum1+sum2+sum3;
}

__device__ fptr_t d_p_ZAKHAROV = d_ZAKHAROV;


/********* QING ************/

__device__ float d_QING(float *x, int n) {
	float sum = 0;
	for(int i = 1; i <= n; i++) {
		sum += powf(powf(x[i-1], 2) - i,2);
	}
	return sum;
}

float h_QING(float *x, int n) {
	float sum = 0;
	for(int i = 1; i <= n; i++) {
		sum += powf(powf(x[i-1], 2) - i,2);
	}
	return sum;
}

__device__ fptr_t d_p_QING = d_QING;

/********* POWELLSINGULAR ************/

__device__ float d_POWELLSINGULAR(float *x, int n) {
	float sum = 0;
	for(int i = 1; i <= n/4; i++) {		
		sum += powf(x[4*i-4] + 10 * x[4*i-3],2) + 
				5 * powf(x[4*i-2] - x[4*i-1],2) +
				powf(x[4*i-3] - x[4*i-2],4) +
				10*powf(x[4*i-4] - x[4*i-1],4);
	}
	return sum;
}

float h_POWELLSINGULAR(float *x, int n) {
	float sum = 0;
	for(int i = 1; i <= n/4; i++) {		
		sum += powf(x[4*i-4] + 10 * x[4*i-3],2) + 
				5 * powf(x[4*i-2] - x[4*i-1],2) +
				powf(x[4*i-3] - x[4*i-2],4) +
				10*powf(x[4*i-4] - x[4*i-1],4);
	}
	return sum;
}

__device__ fptr_t d_p_POWELLSINGULAR = d_POWELLSINGULAR;

/********* POWELLSUM ************/

__device__ float d_POWELLSUM(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(fabs(x[i]), i+2);
	}
	return sum;
}

float h_POWELLSUM(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(fabs(x[i]), i+2);
	}
	return sum;
}

__device__ fptr_t d_p_POWELLSUM = d_POWELLSUM;

/********* MISHRA2 ************/

__device__ float d_MISHRA2(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		sum += x[i] + x[i+1];
	}
	sum = 0.5 * sum;
	return powf(1 + n - sum, n - sum);
}

float h_MISHRA2(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		sum += x[i] + x[i+1];
	}
	sum = 0.5 * sum;

	return powf(1 + n - sum, n - sum);
}

__device__ fptr_t d_p_MISHRA2 = d_MISHRA2;


/********* MISHRA1 ************/

__device__ float d_MISHRA1(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		sum += x[i];
	}
	

	return powf(1 + n - sum, n - sum);
}

float h_MISHRA1(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		sum += x[i];
	}
	

	return powf(1 + n - sum, n - sum);
}

__device__ fptr_t d_p_MISHRA1 = d_MISHRA1;

/********* EXPONENTIAL ************/

__device__ float d_EXPONENTIAL(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(x[i],2);
	}
	sum = -0.5 * sum;

	return -expf(sum);
}

float h_EXPONENTIAL(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(x[i],2);
	}
	sum = -0.5 * sum;

	return -expf(sum);
}

__device__ fptr_t d_p_EXPONENTIAL = d_EXPONENTIAL;

/********* BROWN ************/

__device__ float d_BROWN(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		sum += powf(powf(x[i+1],2), powf(x[i],2) + 1) + powf(powf(x[i],2), powf(x[i+1],2) + 1);
	}	
	return sum;
}

float h_BROWN(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n-1; i++) {
		sum += powf(powf(x[i+1],2), powf(x[i],2) + 1) + powf(powf(x[i],2), powf(x[i+1],2) + 1);
	}	
	return sum;
}

__device__ fptr_t d_p_BROWN = d_BROWN;

/********* ALPINE1 ************/

__device__ float d_ALPINE1(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += fabs(x[i] * sinf(x[i]) + 0.1 * x[i]);
	}	
	return sum;
}

float h_ALPINE1(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += fabs(x[i] * sinf(x[i]) + 0.1 * x[i]);
	}	
	return sum;
}

__device__ fptr_t d_p_ALPINE1 = d_ALPINE1;



/********* CHUNG ************/

__device__ float d_CHUNG(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(x[i], 2);
	}	
	return powf(sum, 2);
}

float h_CHUNG(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(x[i], 2);
	}	
	return powf(sum, 2);
}

__device__ fptr_t d_p_CHUNG = d_CHUNG;


/********* SPHERE ************/

__device__ float d_SPHERE(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(x[i], 2);
	}	
	return sum;
}

float h_SPHERE(float *x, int n) {
	float sum = 0;
	for(int i = 0; i < n; i++) {
		sum += powf(x[i], 2);
	}	
	return sum;
}

__device__ fptr_t d_p_SPHERE = d_SPHERE;


/********* ACKLEY1 ************/

__device__ float d_ACKLEY1(float *x, int n) {
	float sum1 = 0;
	float sum2 = 0;

	for(int i = 0; i < n; i++) {
		sum1 += powf(x[i], 2);
		sum2 += cosf(2*M_PI*x[i]);
	}

	return (-20 * expf(-0.2 * sqrtf(1.0/n * sum1)) - expf(1.0/n * sum2) + 20 + expf(1));
}

float h_ACKLEY1(float *x, int n) {
	float sum1 = 0;
	float sum2 = 0;

	for(int i = 0; i < n; i++) {
		sum1 += powf(x[i], 2);
		sum2 += cosf(2*M_PI*x[i]);
	}

	return (-20 * expf(-0.2 * sqrtf(1.0/n * sum1)) - expf(1.0/n * sum2) + 20 + expf(1));
}

__device__ fptr_t d_p_ACKLEY1 = d_ACKLEY1;


/********* RASTRIGIN ************/

__device__ float d_RASTRIGIN(float *x, int n) {
	float sum = 10 * n;
	for(int i=0; i < n; i++) {
		sum += powf(x[i],2) - 10 * cosf(2 * M_PI * x[i]);
	}

	return sum;
}

float h_RASTRIGIN(float *x, int n) {
	float sum = 10 * n;
	for(int i=0; i < n; i++) {
		sum += powf(x[i],2) - 10 * cosf(2 * M_PI * x[i]);
	}

	return sum;
}

__device__ fptr_t d_p_RASTRIGIN = d_RASTRIGIN;

/********* DIXON ************/

__device__ float d_DIXON(float *x, int n) {
	float sum = powf(x[0] - 1, 2);
	for(int i=1; i < n; i++) {
		sum += (i+1) * powf( 2*powf(x[i],2) - x[i-1], 2);
	}

	return sum;
}

float h_DIXON(float *x, int n) {
	float sum = powf(x[0] - 1, 2);
	for(int i=1; i < n; i++) {
		sum += (i+1) * powf( 2*powf(x[i],2) - x[i-1], 2);
	}

	return sum;
}

__device__ fptr_t d_p_DIXON = d_DIXON;

/********* ROSENBROCK ************/

__device__ float d_ROSENBROCK(float *x, int n) {
	float sum = 0;
 
	for(int i=0; i < n-1; i++) {
		sum += 100*powf(x[i+1] - powf(x[i],2),2) + powf(x[i] - 1, 2);
	}

	return sum;	
}

float h_ROSENBROCK(float *x, int n) {
	float sum = 0;

	for(int i=0; i < n-1; i++) {
		sum += 100*powf(x[i+1] - powf(x[i],2),2) + powf(x[i] - 1, 2);
	}

	return sum;	
}

__device__ fptr_t d_p_ROSENBROCK = d_ROSENBROCK;

/********* STYBLINSK ************/

__device__ float d_STYBLINSK(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n; i++) {
		sum += powf(x[i],4) - 16*powf(x[i],2) + 5*x[i];
	}

	return sum/2;	
}

float h_STYBLINSK(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n; i++) {
		sum += powf(x[i],4) - 16*powf(x[i],2) + 5*x[i];
	}

	return sum/2;	
}

__device__ fptr_t d_p_STYBLINSK = d_STYBLINSK;

/********* SCHWEFEL ************/

__device__ float d_SCHWEFEL(float *x, int n) {
	float sum = 418.9829 * n;
	for(int i=0; i < n; i++) {
		sum += x[i] * sinf(sqrtf(fabs(x[i])));
	}
	if (sum < 0) sum = 0;
	return sum;	
}

float h_SCHWEFEL(float *x, int n) {
	float sum = 418.9829 * n;
	for(int i=0; i < n; i++) {
		sum += x[i] * sinf(sqrtf(fabs(x[i])));
	}
	
	if (sum < 0) sum = 0;
	return sum;	
}

__device__ fptr_t d_p_SCHWEFEL = d_SCHWEFEL;

/********* SALOMON ************/

__device__ float d_SALOMON(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n; i++) {
		sum += powf(x[i],2);
	}
	sum = sqrtf(sum);
	return (1 - cosf(2.0f*M_PI*sum) + 0.1f * sum);
}

float h_SALOMON(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n; i++) {
		sum += powf(x[i],2);
	}
	sum = sqrtf(sum);
	return (1 - cosf(2.0f*M_PI*sum) + 0.1f * sum);
}

__device__ fptr_t d_p_SALOMON = d_SALOMON;

/********* GRIEWANK ************/

__device__ float d_GRIEWANK(float *x, int n) {
	float sum1 = 0;
	for(int i=0; i < n; i++) {
		sum1 += pow(x[i],2);
	}
	sum1 = sum1/4000;

	float sum2 = 1;
	for(int i=0; i < n; i++) {
		sum2 *= cosf(x[i]/sqrtf(i+1));
	}

	return (1 + sum1 - sum2);	
}

float h_GRIEWANK(float *x, int n) {
	float sum1 = 0;
	for(int i=0; i < n; i++) {
		sum1 += pow(x[i],2);
	}
	sum1 = sum1/4000;

	float sum2 = 1;
	for(int i=0; i < n; i++) {
		sum2 *= cosf(x[i]/sqrtf(i+1));
	}

	return (1 + sum1 - sum2);	
}

__device__ fptr_t d_p_GRIEWANK = d_GRIEWANK;

/********* SCHAFER_F6 ************/

__device__ float d_SCHAFER_F6(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n-1; i++) {
		sum += 0.5f+((powf(sinf(sqrtf(powf(x[i], 2) + powf(x[i+1], 2))),2)-0.5f)/powf(1 + 0.001f*(powf(x[i], 2) + 
					powf(x[i+1], 2)), 2));
	}
	return sum;
}

float h_SCHAFER_F6(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n-1; i++) {
		sum += 0.5f+((powf(sinf(sqrtf(powf(x[i], 2) + powf(x[i+1], 2))),2)-0.5f)/powf(1 + 0.001f*(powf(x[i], 2) + 
					powf(x[i+1], 2)), 2));
	}
	return sum;
}

__device__ fptr_t d_p_SCHAFER_F6 = d_SCHAFER_F6;


/********* F25aF28 ************/

__device__ float d_F25aF28(float *x, int n) {
	float sum1 = 0, sum2 = 0;
	for(int i=0; i < n; i++) {
		sum1 += powf(x[i], 2)/powf(2, i);
		if(i<n-1)
			sum2 += powf(x[i+1]*x[i], 2)/powf(2, i+1);
	}
	return sum1+sum2;
}

float h_F25aF28(float *x, int n) {
	float sum1 = 0, sum2 = 0;
	for(int i=0; i < n; i++) {
		sum1 += powf(x[i], 2)/powf(2, i);
		if(i<n-1)
			sum2 += powf(x[i+1]*x[i], 2)/powf(2, i+1);
	}
	return sum1+sum2;
}

__device__ fptr_t d_p_F25aF28 = d_F25aF28;

/********* LEVY ************/

__device__ float d_LEVY(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n-1; i++) {
		sum += powf(x[i]-1,2) * (1+10*powf(sinf(M_PI*x[i]),2));
	}
	return powf(sinf(M_PI * x[0]),2) + sum + powf(x[n-1]-1,2) * (10*powf(sinf(M_PI*x[n-1]),2));
}

float h_LEVY(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n-1; i++) {
		sum += powf(x[i]-1,2) * (1+10*powf(sinf(M_PI*x[i]),2));
	}
	return powf(sinf(M_PI * x[0]),2) + sum + powf(x[n-1]-1,2) * (10*powf(sinf(M_PI*x[n-1]),2));
}

__device__ fptr_t d_p_LEVY = d_LEVY;

/********* PICCIONI ************/

__device__ float d_PICCIONI(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n-1; i++) {
		sum += powf(x[i]-1,2) * (1+10*powf(sinf(M_PI*x[i+1]),2));
	}
	return 10*powf(sinf(M_PI * x[0]),2) + sum + powf(x[n-1]-1,2);
}

float h_PICCIONI(float *x, int n) {
	float sum = 0;
	for(int i=0; i < n-1; i++) {
		sum += powf(x[i]-1,2) * (1+10*powf(sinf(M_PI*x[i+1]),2));
	}
	return 10*powf(sinf(M_PI * x[0]),2) + sum + powf(x[n-1]-1,2);
}

__device__ fptr_t d_p_PICCIONI = d_PICCIONI;

/********* CLUSTER MARIO *********/ // Como acessar points dessa funcao [BRUNO]
__device__ float* g_d_points; //esse é o global device
//__constant__ int g_d_dim; 
	//__constant__ int g_d_npoints; 

float* points;
const float my_pi = 3.14;


__device__ float d_CLUSTER(float *x, int n) {
	float sum = 0;
	return g_d_points[0] ; //[BRUNO] ERRO está aqui!
}

float h_CLUSTER(float *x, int n) {
	float sum = 0;

	return points[0];
}

__device__ fptr_t d_p_CLUSTER = d_CLUSTER;



int main(int argc, char **argv){
	fptr_t d_obj_f, h_obj_f;
	int n = 25;
	double cutoff_time = 60;
	float hs = 0.5, he = 0.0001, ep = 1/powf(2,13);
	int gpu = 1;
	int max_points = 128;
	int seed = time(NULL);

	Results res1, res2, res3, res4;



	int n_points;
	int dim;


	std::string s;
	
	std::getline(std::cin, s);
	std::stringstream st(s);
	st >> n_points >> dim;
	points = new float[n_points*dim];
	n = dim;

	float CLUSTER_L[n], CLUSTER_U[n];

	for(int i=0;i<n;++i) {
		CLUSTER_L[i] = 1000;
		CLUSTER_U[i] = -1000;
	}

	for(int i=0;i<n_points;++i) {
		
		std::getline(std::cin, s);
		std::stringstream st(s);
		for(int d=0; d<dim; ++d) {
			st >> points[4*i + d];

			if(CLUSTER_L[d] > points[4*i + d])
				CLUSTER_L[d] = points[4*i + d]; //Min da dim
			if(CLUSTER_U[d] < points[4*i + d])
				CLUSTER_U[d] = points[4*i + d]; //Max da dim

		}



	}
        float *d_points; // 


        cudaMalloc( &d_points, n_points * dim * sizeof(float)); //aloco e copio
        cudaMemcpy(d_points, points, n_points * dim * sizeof(float), cudaMemcpyHostToDevice);

        cudaMemcpyToSymbol(g_d_points, &d_points, sizeof(d_points));



       //   cudaMemcpyToSymbol(g_d_dim, &dim, sizeof(int));
       // cudaMemcpyToSymbol(g_d_npoints, &n_points, sizeof(int));

	/*for(int i=0;i<n_points;i++){
		for(int d=0;d<dim;d++)
			std::cout << points[i][d] << " ";
		std::cout << std::endl;
	}*/

        //gpuErrchk(cudaMalloc(&d_points, n_points* dim * sizeof(float)));



	
	/*************** CLUSTER ***********************/


	cudaMemcpyFromSymbol(&d_obj_f, d_p_CLUSTER, sizeof (fptr_t));
	h_obj_f = h_CLUSTER;

	// // float POWELLSINGULAR_L[n], POWELLSINGULAR_U[n];
	// // for (int i = 0; i < n; i++) {
	// // 	POWELLSINGULAR_L[i] = -4;
	// // 	POWELLSINGULAR_U[i] = 5;
	// // }


	//res1 = gpu_cgrasp(n, CLUSTER_L, CLUSTER_U, points, n_points, dim, 0, seed, ep, hs, he, cutoff_time,
	 //			   	  max_points, 0, h_obj_f, d_obj_f);

	printf("total time: %lf\n", res1.time);
	res2 = gpu_cgrasp(n, CLUSTER_L, CLUSTER_U, points, n_points, dim, 0.0, seed, ep, hs, he, cutoff_time,
				   	  max_points, 1, h_obj_f, d_obj_f);

	//res3 = pcgrasp(n, POWELLSINGULAR_L, POWELLSINGULAR_U, points, n_points, dim, 0.0, seed, ep, hs, he, cutoff_time,
	 //			   	  max_points, 0, h_obj_f, d_obj_f);

	// res4 = pcgrasp(n, POWELLSINGULAR_L, POWELLSINGULAR_U, 0.0, seed, ep, hs, he, cutoff_time,
	//  			   	  max_points, 1, h_obj_f, d_obj_f);

	printf("\n\nInstance: CLUSTER\n");
	printf("total time: %lf\n", res1.time);
	printf("total time (gpu): %lf\n", res2.time);
	printf("total time (par): %lf\n", res3.time);
	printf("total time (gpar): %lf\n", res4.time);
	printf("total evaluations: %ld\n", res1.evaluations);
	printf("total evaluations (gpu): %ld\n", res2.evaluations);
	printf("total evaluations (par): %ld\n", res3.evaluations);
	printf("total evaluations (gpar): %ld\n", res4.evaluations);
	printf("error: %lf\n", fabs(res1.best));
	printf("error (gpu): %lf\n", fabs(res2.best));
	printf("error (par): %lf\n", fabs(res3.best));
	printf("error (gpar): %lf\n", fabs(res4.best));



	return 0;
}
