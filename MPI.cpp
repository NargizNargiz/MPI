#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
using namespace std;

int main(int argc, char **argv){
        int size, MyP, N, M;
        double  temp;
        double t0,t1;
        MPI_Status stat;
        MPI_Status stat1;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &MyP);
 //     cout << argv[1] << endl;
 //     ifstream file(argv[1]);
 //     cout << argv[1] <<endl;
 //     file >> M;
        M = 256;
        cout << M << endl;
        N = M/size;// size = 4, N = 2, M = 4
        double **matrix = new double *[N];
        double **E = new double *[N];
        double *V = new double [M];
        double *V1 = new double [M];
        for (int i = 0; i < N; i++)
                matrix[i] = new double [M];
        for (int i = 0; i < N; i++)
                E[i] = new double [M];
        srand(time(NULL));
        for(int i = 0; i < N; i++)
                for(int j = 0; j < M; j++)
                        matrix[i][j] =  rand() % 101;
        for(int i = 0; i < N; i++){
                for(int j = 0; j < M; j++){
                 if (N*MyP + i == j)
                        E[i][j] = 1.0;
                else
                        E[i][j] = 0.0;
                }
        }
        t0 = MPI_Wtime();
        for(int p = 0; p < size; p++){
                for(int k = 0; k < N; k++){
                        if(MyP == p){
                                temp = 1.0/matrix[k][N*p +k];
                                for (int i = 0; i < M; i++){
                                        matrix[k][i] = matrix[k][i]*temp;
                                        E[k][i] = E[k][i]*temp;
                                }
                                //отправляем всем остальнымж
                                for (int m = p+1; m < size; m++){
                                         MPI_Send(&matrix[k][0], M, MPI_DOUBLE, m, 1,MPI_COMM_WORLD);
                                         MPI_Send(&E[k][0], M, MPI_DOUBLE, m, 1,MPI_COMM_WORLD);
                                }//обнуляем оставшиеся строки в подматрице данного процессора
                                for (int i = k+1; i < N; i++){
                                        for (int j = 0 ; j< M; j++){
                                                matrix[i][j] = matrix[i][j]-matrix[i][N*p +k] * matrix[k][j];
                                                E[i][j] = E[i][j] - E[i][N*p +k] * matrix[k][j];
                                        }
                                }
                        }else{
                                if (MyP > p){
                                        MPI_Recv(V, M, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &stat);
                                        MPI_Recv(V1, M, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &stat1);
                                        for (int i = 0; i < N; i++)
                                                for (int j = 0;j < M; j++){
                                                        matrix[i][j] = matrix[i][j] - matrix[i][N*p+k]*V[j];
                                                        E[i][j] = E[i][j] - E[i][N*p +k]*V[j];
                                                }
                                }
                        }
                }
        }                                                                                                                                                     
        for( int p = size-1; p >= 0; p--){
             for(int k = N-1; k >= 0; k--){
                        if(MyP == p){
                                for (int m = p-1 ; m >= 0; m--){
                                        MPI_Send(&matrix[k][0], M, MPI_DOUBLE, m, 1, MPI_COMM_WORLD);
                                        MPI_Send(&E[k][0], M, MPI_DOUBLE, m, 1, MPI_COMM_WORLD);
                                }
                                for (int i = k-1; i>= 0; i--)
                                        for (int j = 0; j < M; j++ ){
                                                matrix[i][j] = matrix[i][j] - matrix[i][N*p+k] *matrix[k][j];
                                                E[i][j] = E[i][j] - E[i][N*p+k] *matrix[k][j];
                                        }
                        }else{
                                if (MyP < p){
                                        MPI_Recv(V, M, MPI_DOUBLE, p, 1,MPI_COMM_WORLD, &stat);
                                        MPI_Recv(V1, M, MPI_DOUBLE, p, 1,MPI_COMM_WORLD, &stat1);
                                        for (int i = N - 1; i >= 0; i--)
                                                for (int j = 0; j < M; j++){
                                                        matrix[i][j]= matrix[i][j] - matrix[i][N*p+k]*V[j];
                                                        E[i][j]= E[i][j] - E[i][N*p+k]*V[j];
                                                }
                                }
                        }
                }
        }   
        t1 = MPI_Wtime()-t0;
        cout << t1 <<endl;

        for (int i = 0; i < N; i++){
                delete [] E[i];
        }
        delete [] E;
        for (int i = 0; i < N; i++)
                delete [] matrix[i];
        delete [] matrix;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return(0);
}
