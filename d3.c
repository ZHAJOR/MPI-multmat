#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
//#include <omp.h>

void fillMatrice(int*** matrice, int matrice_size){
	*matrice = malloc(matrice_size * sizeof(int*));
	for(int i=0; i<matrice_size; ++i){
		(*matrice)[i] = malloc(matrice_size * sizeof(int));
		for(int j=0; j<matrice_size; ++j){
			scanf("%d", &(*matrice)[i][j]);
		}
	}
}
int** initMatrice(int matrice_size){
	int** matrice = malloc(matrice_size * sizeof(int*));
	for(int i=0; i<matrice_size; ++i){
		matrice[i] = malloc(matrice_size * sizeof(int));
		for(int j=0; j<matrice_size; ++j){
			matrice[i][j] = 0;
		}
	}
	return matrice;
}

void printMatrice(int** matrice, int matrice_size){
	for(int i=0; i<matrice_size; ++i){
		for(int j=0; j<matrice_size; ++j){
			printf("%5d ", matrice[i][j]);
		}
		printf("\n");
	}
}


void arrayBroadcast(int** a, int matrice_size, int master){
	for(int i = 0; i < matrice_size; ++i){
		MPI_Bcast(a[i], matrice_size, MPI_INT, master, MPI_COMM_WORLD);
	}
}

int main (int argc, char *argv[])
{
	int id, p;
	int matrice_size, from, to;
	int** a;
	int** b;
	double elapsed_time;
	int *rcvcounts;
	int *displs;
	int *positions;

	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	rcvcounts = calloc(p, sizeof(int));
	displs = calloc(p, sizeof(int));
	positions = calloc(p, sizeof(int));

       	if(id == 0){
		scanf("%d", &matrice_size);
		fillMatrice(&a, matrice_size);
		fillMatrice(&b, matrice_size);
		elapsed_time = -MPI_Wtime();

		p = p>matrice_size ? matrice_size : p;
		int modulo = matrice_size % p;
		int nb = 0;
		for(int i = 0; i < p ;++i){
			int numberElements = matrice_size/p;
			if(modulo > 0){
				++numberElements;
				--modulo;
			}
			if(i == 0)
				positions[i] = numberElements;
			else
				positions[i] = positions[i-1] + numberElements;
			numberElements *= matrice_size;
			rcvcounts[i] = numberElements;
			displs[i] = nb;
			nb += numberElements;
		}
	}
	MPI_Bcast(&matrice_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(rcvcounts, p, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, p, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(positions, p, MPI_INT, 0, MPI_COMM_WORLD);
	if(id != 0){
		a=initMatrice(matrice_size);
		b=initMatrice(matrice_size);
	}
	int *c = malloc(sizeof(int) * matrice_size * matrice_size);
	arrayBroadcast(a, matrice_size, 0);
	arrayBroadcast(b, matrice_size, 0);

	to = positions[id];
	from = to - (rcvcounts[id] / matrice_size);

//	#pragma omp parallel for
	for(int i =from; i<to; ++i) {
		for(int j = 0; j<matrice_size; ++j){
			int su = 0;
			for(int k=0; k<matrice_size; ++k){
				su += a[i][k] * b[k][j];
			}
			c[i *matrice_size + j] = su;
		}
	}

	MPI_Gatherv(&c[from * matrice_size], rcvcounts[id], MPI_INT, c, rcvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);


	if(id == 0 && argc > 1 && argv[1][0] == '-' && argv[1][1] == 'p'){
		printf("-------------\n");
		for(int i=0; i<matrice_size; ++i){
			for(int j=0; j<matrice_size; ++j){
				printf("%5d ", c[i * matrice_size + j]);
			}
			printf("\n");
		}
		printf("-------------\n");
	}
	elapsed_time +=MPI_Wtime();
	//if(!id){
		//printf("%10.6f\n", elapsed_time);
	//}

	MPI_Finalize();
	return 0;
}
