#include "page_rank.h"

#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"


// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double* solution, double damping, double convergence)
{

    // initialize vertex weights to uniform probability. Double
    // precision scores are used to avoid underflow for large graphs

    int numNodes = num_nodes(g);
    double equal_prob = 1.0 / numNodes;
    double sum_sinks = 0.0;
    double* score_new = new double[numNodes];
    double * score_old = new double[numNodes];
    #pragma omp parallel for reduction(+:sum_sinks) if(omp_get_max_threads() > 1)
    for (int i = 0; i < numNodes; ++i) {
        score_old[i] = equal_prob;
        if (outgoing_begin(g, i) == outgoing_end(g, i)) {
            sum_sinks += score_old[i];
        }
    }

    bool converged = false;
    while (!converged) {
        double global_diff = 0.0;
        double new_sum_sinks = 0.0;

        // compute score_new[vi] for all nodes vi:
        #pragma omp parallel for reduction(+:global_diff, new_sum_sinks) schedule(dynamic, 16) if(omp_get_max_threads() > 1)
        for (int i = 0; i < numNodes; i++) {
            score_new[i] = 0.0;
            // iterate through all incoming edges of node i
            const Vertex* in_begin = incoming_begin(g, i);
            const Vertex* in_end = incoming_end(g, i);
            for (const Vertex* v = in_begin; v != in_end; ++v) {
                score_new[i] += score_old[*v] / outgoing_size(g, *v);
            }
            score_new[i] = (damping * score_new[i]) + (1.0 - damping) /(numNodes);
            score_new[i] += damping * sum_sinks / (numNodes);
            if (outgoing_size(g, i) == 0) {
                new_sum_sinks += score_new[i];
            }
            // compute how much per-node scores have changed
            global_diff += fabs(score_new[i] - score_old[i]);
        }

        // update solution
        std::swap(score_new, score_old);

        sum_sinks = new_sum_sinks;

        // quit once algorithm has converged
        converged = (global_diff < convergence);
    }
    #pragma omp parallel for reduction(+:sum_sinks) if(omp_get_max_threads() > 1)
    for (int i = 0; i < numNodes; ++i) {
        solution[i] = score_old[i];
    }

    /*
        CS149 students: Implement the page rank algorithm here.  You
        are expected to parallelize the algorithm using openMP.  Your
        solution may need to allocate (and free) temporary arrays.

        Basic page rank pseudocode is provided below to get you started:

        // initialization: see example code above
        score_old[vi] = 1/numNodes;

        while (!converged) {

        // compute score_new[vi] for all nodes vi:
        score_new[vi] = sum over all nodes vj reachable from incoming edges
                            { score_old[vj] / number of edges leaving vj  }
        score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

        score_new[vi] += sum over all nodes v in graph with no outgoing edges
                            { damping * score_old[v] / numNodes }

        // compute how much per-node scores have changed
        // quit once algorithm has converged

        global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
        converged = (global_diff < convergence)
        }

    */
}
