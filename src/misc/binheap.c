/* binheap.c */

/* S. Engblom 2017-02-16 */
/* J. Cullhed 2008-08-04. */

/* Binary heap needed by the NRM and NSM. */

/* This is an indexed heap. If node is the index of an element in
   "data", the place in the heap is given by INDEX[node]. The place in
   the "data" held by node n in the heap is given by INDEX2[INDEX[n]].
*/

#include <R_ext/Visibility.h>

static void
percolate_down(int n1, double *data, int *INDEX, int *INDEX2, int N)
{
    int child;
    int node = n1;
    const double key = data[node];
    const int j = INDEX[node];

    while ((child = (node << 1) + 1) < N) {
        if (child != N - 1 && data[child + 1] < data[child])
            child++;

        if (data[child] < key) {
            data[node] = data[child];
            const int rxn = INDEX[child];
            INDEX2[rxn] = node;
            INDEX[node] = rxn;
        } else {
            break;
        }

        node = child;
    }

    data[node] = key;
    INDEX[node] = j;
    INDEX2[j] = node;
}

static void percolate_up(int node, double *data, int *INDEX, int *INDEX2)
{
    int parent;
    const double key = data[node];
    const int j = INDEX[node];

    do {
        parent = (node - 1) >> 1;

        if (key < data[parent]) {
            data[node] = data[parent];
            const int rxn = INDEX[parent];
            INDEX2[rxn] = node;
            INDEX[node] = rxn;
        } else {
            break;
        }

        node = parent;
    } while (parent > 0);

    data[node] = key;
    INDEX[node] = j;
    INDEX2[j] = node;
}

attribute_hidden
    void initialize_heap(double *data, int *INDEX, int *INDEX2, int N)
{
    for (int i = (N - 1) >> 1; i >= 0; i--)
        percolate_down(i, data, INDEX, INDEX2, N);
}

attribute_hidden
    void update(int node, double *data, int *INDEX, int *INDEX2, int N)
{
    const int parent = (node - 1) >> 1;

    if (node > 0 && data[node] < data[parent])
        percolate_up(node, data, INDEX, INDEX2);
    else
        percolate_down(node, data, INDEX, INDEX2, N);
}
