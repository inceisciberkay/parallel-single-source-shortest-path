#include "util.h"
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void usage() {
  printf("Usage: mpirun -np <number_of_processors> ./sssp_parallel "
         "<input_file> <source_vertex> <output_file>");
}

int min(int a, int b) { return a > b ? b : a; }

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Incorrect number of arguments.\n");
    usage();
    return 1;
  }

  char *infile_name = argv[1];
  int source_vertex = atoi(argv[2]);
  char *outfile_name = argv[3];

  int myid, numprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int num_vertices, num_edges, *offsets, *edges, *weights;

  if (myid == 0) { // Master Process
    // reading input file and initializing arrays
    int success = read_file(infile_name, &num_vertices, &num_edges, &offsets);

    if (!success) {
      fprintf(stderr, "Input file could not be read");
      return 1;
    }

    edges = &(offsets[num_vertices + 1]);
    weights = &(offsets[num_vertices + 1 + num_edges]);

    // sending num_vertices information to worker processes
    // since only the master process reads the input file, only it has
    // num_vertices information. Therefore, this information must be shared with
    // the worker processes.
    for (int worker_id = 1; worker_id < numprocs; worker_id++) {
      MPI_Send(&num_vertices, 1, MPI_INT, worker_id, 0, MPI_COMM_WORLD);
    }
  } else { // Worker Processes
    // receiving num_vertices information from master process
    MPI_Recv(&num_vertices, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double start = MPI_Wtime();

  // partitioning vertices
  int *send_counts = malloc(numprocs * sizeof(int));
  int *displs = malloc(numprocs * sizeof(int));

  // actually, there is no need to perform this check (myid == 0), but I keep it
  // to eliminate unnecessary computations in the worker processes
  if (myid == 0) {
    // setting send_counts and displs arrays
    for (int i = 0; i < numprocs; i++) {
      send_counts[i] = num_vertices / numprocs + 1;
      displs[i] = i * num_vertices / numprocs;
    }
  }

  int rcvsize = num_vertices / numprocs + 1;
  int *local_vertices = malloc(rcvsize * sizeof(int));
  MPI_Scatterv(offsets, send_counts, displs, MPI_INT, local_vertices, rcvsize,
               MPI_INT, 0, MPI_COMM_WORLD);

  // partitioning edges and weights
  if (myid == 0) { // only root process has 'offsets' array initialized from
                   // the input file
    // setting send_counts and displs arrays
    for (int i = 0; i < numprocs; i++) {
      send_counts[i] = offsets[(i + 1) * (num_vertices / numprocs)] -
                       offsets[i * (num_vertices / numprocs)];
      if (i == 0) {
        displs[i] = 0;
      } else {
        displs[i] = displs[i - 1] + send_counts[i - 1];
      }
    }
  }

  rcvsize = local_vertices[rcvsize - 1] - local_vertices[0];
  int *local_edges = malloc(rcvsize * sizeof(int));
  int *local_weights = malloc(rcvsize * sizeof(int));
  MPI_Scatterv(edges, send_counts, displs, MPI_INT, local_edges, rcvsize,
               MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(weights, send_counts, displs, MPI_INT, local_weights, rcvsize,
               MPI_INT, 0, MPI_COMM_WORLD);

  // performing computation
  int *D = malloc(num_vertices * sizeof(int));
  for (int i = 0; i < num_vertices; i++) {
    D[i] = INT32_MAX;
  }
  D[source_vertex] = 0;

  int *R = malloc((num_vertices / numprocs) * sizeof(int));

  int changed = 1; // for early termination
  for (int j = 0; j < num_vertices - 1 && changed; j++) {
    changed = 0;
    // displacement of local_edges and local_weigths
    int displ = local_vertices[0];
    for (int i = 0; i < num_vertices / numprocs; i++) {
      int single_min = INT32_MAX;
      for (int k = local_vertices[i]; k < local_vertices[i + 1]; k++) {
        single_min = min(single_min, D[local_edges[k - displ]] == INT32_MAX
                                         ? INT32_MAX // prevent overflow
                                         : local_weights[k - displ] +
                                               D[local_edges[k - displ]]);
      }
      R[i] = min(D[i + myid * (num_vertices / numprocs)], single_min);
      if (R[i] != D[i + myid * (num_vertices / numprocs)]) {
        changed = 1;
      }
    }
    MPI_Allgather(R, num_vertices / numprocs, MPI_INT, D,
                  num_vertices / numprocs, MPI_INT, MPI_COMM_WORLD);

    // check if all processes have "changed = 1":
    int sum = 0;
    MPI_Allreduce(&changed, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    changed = sum == 0 ? 0 : 1;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double end = MPI_Wtime();
  double time = end - start;
  printf("Time: %f seconds\n", time);

  // outputting the result
  if (myid == 0) {
    printResults(outfile_name, D, num_vertices);
  }

  free(R);
  free(D);
  free(local_weights);
  free(local_edges);
  free(local_vertices);
  free(displs);
  free(send_counts);

  MPI_Finalize();

  return 0;
}
