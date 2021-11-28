#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <iostream>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    if (frontier->count == 0) {
        return;
    }
    int next_dist = distances[frontier->vertices[0]] + 1;
    #pragma omp parallel for schedule(dynamic, 128) if(omp_get_max_threads() > 1)
    for (int i=0; i<frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            bool updated = __sync_bool_compare_and_swap(
                &distances[outgoing],
                NOT_VISITED_MARKER,
                next_dist);
            
            if (updated) {
                int index;
                #pragma omp atomic capture
                index = new_frontier->count++;
                new_frontier->vertices[index] = outgoing;
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

int bottom_up_step(
    Graph g,
    int* distances,
    int curr_dist)
{
    int new_frontier_size = 0;

    #pragma omp parallel for reduction(+:new_frontier_size) schedule(dynamic, 128)
    for (int i = 0; i < g->num_nodes; i++) {
        // for each unvisited vertex, iterate through the incoming edges
        // to see if any of them are in the frontier (i.e. sol
        // dist is curr_dist; if so, add them to the new frontier
        // (i.e. set sol dist to curr_dist + 1).
        if (distances[i] != NOT_VISITED_MARKER)
            continue;
        const Vertex* inc_beg = incoming_begin(g, i);
        const Vertex* inc_end = incoming_end(g, i);
        for (const Vertex* u = inc_beg; u != inc_end; u++) {
            if (distances[*u] == curr_dist) {
                new_frontier_size++;
                distances[i] = curr_dist + 1;
                break;
            }
        }
    }
    return new_frontier_size;
}


void bfs_bottom_up(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    for (int i = 0; i < graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    int curr_dist = 0;
    sol->distances[ROOT_NODE_ID] = 0;

    while (true) {
        int new_frontier_size = bottom_up_step(graph, sol->distances, curr_dist);
        if (new_frontier_size == 0) break;
        curr_dist++;
    }
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.

    // Top-down works better when frontier is a small fraction of the
    // total number of nodes, and bottom-up works better when it is
    // a larger fraction. So we can just monitor the ratio of frontier
    // size to total number, and then when it is sufficiently high
    // we use a bottom-up step or when it is low we use a top-down step.

    threshold = 0.33;
    for (int i = 0; i < graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }
    int curr_dist = 0;
    sol->distances[ROOT_NODE_ID] = 0;
    

    while(true) {
        if ((float)frontier_size / (float)graph->num_nodes > threshold) {
            
        }
    }
}
