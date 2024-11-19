#include "util.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void usage() {
  printf("Usage: ./sssp_serial <input_file> <source_vertex> <output_file>");
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

  int num_vertices, num_edges, *offsets, *edges, *weights, *vertices;
  int success = read_file(infile_name, &num_vertices, &num_edges, &offsets);

  if (!success) {
    fprintf(stderr, "Input file could not be read");
    return 1;
  }

  edges = &(offsets[num_vertices + 1]);
  weights = &(offsets[num_vertices + 1 + num_edges]);
  vertices = offsets;

  int *D = malloc(num_vertices * sizeof(int));
  for (int i = 0; i < num_vertices; i++) {
    D[i] = INT32_MAX;
  }
  D[source_vertex] = 0;

  int *R = malloc(num_vertices * sizeof(int));

  clock_t start = clock();
  // computation
  int changed = 1; // for early termination
  for (int j = 0; j < num_vertices - 1 && changed; j++) { // repeat V-1 times
    changed = 0;
    for (int i = 0; i < num_vertices; i++) {
      int single_min = INT32_MAX;
      for (int k = vertices[i]; k < vertices[i + 1]; k++) {
        single_min = min(single_min, D[edges[k]] == INT32_MAX
                                         ? INT32_MAX // prevent overflow
                                         : weights[k] + D[edges[k]]);
      }
      R[i] = min(D[i], single_min);
      if (R[i] != D[i])
        changed = 1;
    }
    // swap two pointers of arrays instead of copying R into D (D <- R)
    int *temp;
    temp = D;
    D = R;
    R = temp;
  }

  clock_t end = clock();
  double time = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Time: %f seconds\n", time);

  printResults(outfile_name, D, num_vertices);

  free(R);
  free(D);
  free(offsets);

  return 0;
}
