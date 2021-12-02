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

using namespace std;

void vertex_set_clear(vertex_set* list) {
    list->total_count = 0;
    // set each count to zero
    for (int i = 0; i < omp_get_max_threads(); i++) {
        list->counts[i] = 0;
    }
}

void vertex_set_init(vertex_set* list, int count) {
    // malloc vertex array pointer list (one for each thread)
    list->vertices = (int**) malloc(sizeof(int*) * omp_get_max_threads());
    // malloc vertex array for each thread
    for (int i = 0; i < omp_get_max_threads(); i++) {
        list->vertices[i] = (int*) malloc(sizeof(int) * count);
    }
    // malloc counts array
    list->counts = (int*) malloc(sizeof(int) * omp_get_max_threads());
    vertex_set_clear(list);
}

void print_vertex_set(vertex_set* list) {
    printf("vertex set: ");
    // for each thread
    for (int i = 0; i < omp_get_max_threads(); i++) {
        printf("[ ");
        // for each vertex in the thread list
        for (int j = 0; j < list->counts[i]; j++) {
            printf("%d ", list->vertices[i][j]);
        }
        printf("], ");
    }
    printf("\n");
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances,
    int curr_dist)
{
    if (frontier->total_count == 0) {
        return;
    }
    int next_dist = curr_dist + 1;
    #pragma omp parallel for if(omp_get_max_threads() > 1) schedule(dynamic)
    for (int i=0; i<frontier->total_count; i++) {

        // get slice
        int slice = 0;
        int remainder = i;
        while (remainder >= frontier->counts[slice]) {
            remainder -= frontier->counts[slice];
            slice++;
        }

        int node = frontier->vertices[slice][remainder];
        int thread_num = omp_get_thread_num();

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            // bool updated = __sync_bool_compare_and_swap(
            //     &distances[outgoing],
            //     NOT_VISITED_MARKER,
            //     next_dist);
            
            // if (updated) {
            if (distances[outgoing] == NOT_VISITED_MARKER) {
                distances[outgoing] = next_dist;
                int index = new_frontier->counts[thread_num]++;
                new_frontier->vertices[thread_num][index] = outgoing;
            }
        }
    }

    for (int j=0; j<omp_get_max_threads(); j++) {
        new_frontier->total_count += new_frontier->counts[j];
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
    frontier->vertices[0][frontier->counts[0]++] = ROOT_NODE_ID;
    frontier->total_count++;
    sol->distances[ROOT_NODE_ID] = 0;
    int curr_dist = 0;

    while (frontier->total_count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances, curr_dist);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;

        curr_dist++;
    }
}

int bottom_up_step(
    Graph g,
    int* distances,
    int curr_dist)
{
    int new_frontier_size = 0;

    #pragma omp parallel for reduction(+:new_frontier_size) schedule(dynamic)
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
    
    float threshold = 0.33;

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
    frontier->vertices[0][frontier->counts[0]++] = ROOT_NODE_ID;
    frontier->total_count++;
    sol->distances[ROOT_NODE_ID] = 0;
    int curr_dist = 0;
    int frontier_size = 1;
    bool using_top_down = true;

    while (frontier_size != 0) {
        if ((float)frontier_size / graph->num_nodes < threshold && !using_top_down) {
            // need to switch to using top down, i.e. need to rebuild the frontier
            using_top_down = true;
            vertex_set_clear(frontier);
            // add all nodes with distance curr_dist to the frontier
            #pragma omp parallel for
            for (int i=0; i<graph->num_nodes; i++) {
                if (sol->distances[i] == curr_dist) {
                    int thread_num = omp_get_thread_num();
                    frontier->vertices[thread_num][frontier->counts[thread_num]++] = i;
                }
            }
            for (int j=0; j<omp_get_max_threads(); j++) {
                frontier->total_count += frontier->counts[j];
            }
        } else if ((float)frontier_size / graph->num_nodes >= threshold) {
            using_top_down = false;
        }

        if (using_top_down) {
            vertex_set_clear(new_frontier);
            top_down_step(graph, frontier, new_frontier, sol->distances, curr_dist);
            // swap pointers
            vertex_set* tmp = frontier;
            frontier = new_frontier;
            new_frontier = tmp;
            frontier_size = frontier->total_count;
        } else {
            frontier_size = bottom_up_step(graph, sol->distances, curr_dist);
        }
        curr_dist++;
    }
}
