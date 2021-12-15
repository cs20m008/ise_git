#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

/***
 * Return seconds elapsed since t0, with t0 = 0 the epoch
 ***/
double stop_watch(double t0)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec/1e6 - t0;
}

/*
 * Function reads `n' lines from file with name `fname[]` and returns
 * each line as a double precision number in `arr[]'. `arr[]' should
 * have enough space for `n' double precision numbers
 */
/*void readlines(int n, double *arr, char fname[])
{
  FILE *fp = fopen(fname, "r");
  double line_val;
  for(int i=0; i<n; i++) {
    fscanf(fp, "%lf", &line_val);
    arr[i] = line_val;
  }
  fclose(fp);
  return;
}*/

int main(int argc, char *argv[])
{
  char hname[256];
  gethostname(hname, 256);
    
  MPI_Init(&argc, &argv);
  
  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*
   * Abort if the number of processes does not divide `nelems' exactly
   */
  int nelems = 360360;
  if(nelems % nproc != 0) {
    fprintf(stderr, " nelems = %d not divisible by nproc = %d\n", nelems, nproc);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  vector<double>x(nelems);
  vector<double>y(nelems);
 /* double *x = malloc(nelems*sizeof(double));
  double *y = malloc(nelems*sizeof(double));*/
  if(rank == 0) {
   for(int i=0;i<nelems;i++)
   {
       x[i]=rand%100;
       y[i]=rand%100;
   }
  }

  /*
   * `nelems_loc': local number of elements, i.e. number of elements
   * pre process
   *
   * `x_loc[]': local `x[]' array
   * `y_loc[]': local `y[]' array
   */
  int nelems_loc = nelems/nproc;
 /* double *x_loc = malloc(nelems_loc*sizeof(double));
  double *y_loc = malloc(nelems_loc*sizeof(double));*/
  vector<double>x_loc(nelems_loc);
  vector<double>y_loc(nelems_loc);

  MPI_Scatter(x, nelems_loc, MPI_DOUBLE, x_loc, nelems_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(y, nelems_loc, MPI_DOUBLE, y_loc, nelems_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double t = stop_watch(0);

  double dot_loc = 0;
#pragma omp parallel for reduction(+:dot_loc)
  for(int i=0; i<nelems_loc; i++)
    dot_loc += x_loc[i]*y_loc[i];

  double dot;
  MPI_Reduce(&dot_loc, &dot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  t = stop_watch(t);

#pragma omp parallel
  {
    int nthr = omp_get_num_threads();
#pragma omp single
    {
      if(rank == 0)
	printf(" Used %d MPI processes, %d OpenMP threads, dot(x,y) = %lf   took %12.8f msec\n", nproc, nthr, dot, t*1e3);
    }
  }

  MPI_Finalize();
  return 0;
}