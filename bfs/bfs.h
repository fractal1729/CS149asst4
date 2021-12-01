#ifndef __BFS_H__
#define __BFS_H__

//#define DEBUG

#include "common/graph.h"

struct solution
{
  int *distances;
};

struct vertex_set {
  // # of vertices in the set
  int total_count;
  // # of vertices in each thread's set
  int *counts;
  // max size of buffer vertices 
  int max_vertices_per_thread;
  // array of arrays of vertex ids
  int **vertices;
};

void bfs_top_down(Graph graph, solution* sol);
void bfs_bottom_up(Graph graph, solution* sol);
void bfs_hybrid(Graph graph, solution* sol);

#endif
