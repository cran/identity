/* nodehash.h */

/* Functions needed to use a hash to store and recall computed probabilities
   for sets of equivalence classes.
*/

#define NID 4
struct hash_cell {
   Probvec_t prob;
   //int connectarr[NGENES];
   int nodeid[NID];
   int cstate;
   int isempty;
   struct hash_cell *next;
   struct hash_cell *prev; /* In the hash table, this points to the last
                              cell in the collision list. */
};
typedef struct hash_cell hash_cell_t;

/**
 * Initialize the hash table.
 *
 * @param size the maximum amount of memory for hash cells that will be
 * allocated during execution of the program.  Given in megabytes.
 */
void hashinit (int memsize);
void hashstore (Ibdgraph_t *node, Probvec_t *prob);
int hashfind (Ibdgraph_t *node, Probvec_t *prob);

void hash_free (void);
