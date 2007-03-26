/* ibdgraph.h */
#ifndef IBDGRAPH_H
#define IBDGRAPH_H

/* Functions for using ibdgraphs. */

#define NGENES    4       /* For generalized kinship coefs up to four people. */
#define NCOEF     15      /* number of detailed identity states for 4 genes.*/
#define NIDSTATE  9       /* number of condensed identity states */

typedef struct  {
   /* all elements are multiplied by 1/2^coef to get prob of that state.*/
    double istate[NCOEF];
} Probvec_t;

typedef struct ibdgraph {
    int genelist[NGENES];        /* ID's of people in the list. */
    int connectarr[NGENES];      /* Identical entries are connected. */
    int connstate; /* Which of the 15 possible connection states connectarr is in. */
} Ibdgraph_t;

/**
 * Initialize ibd graph.
 *
 * @param pnode pointer to Ibdgraph_t
 * @param id1
 * @param id2
 * @param id3
 * @param id4
 */
void
ibdgr_init (Ibdgraph_t *pnode, int id1, int id2, int id3, int id4);

/**
 * The recursive algorithm for following branches up the pedigree for a given
 * list of genes (ie a node), where a gene is a randomly chosen gene from an
 * individual. The individuals are what are listed in geneset.
 *
 * @param pprob
 * @param genenode
 */
void
nodeprob (Probvec_t *pprob, Ibdgraph_t *genenode);

#endif  /* IBDGRAPH_H */
