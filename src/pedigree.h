/* pedigree.h - Version 5 */
#ifndef PEDIGREE_H
#define PEDIGREE_H

typedef struct person_tag {
   int id;                /* Sequential id used internally by the programs. */
   int findiv;    /* The id defined externally by the user for this person. */
   int parent[2];    // parent[0] id number of mother, parent[1] id of father.
} Person_t;

struct idmap {
   int id;                      // internal id.
   int fid;                     // external id (id used in input file).
};
typedef struct idmap idmap_t;

/* General header file for pedigree structures. */

void create_pedigree (int n, int *fid, int *ma, int *pa);

int getparent(int whichpar, int indiv);

int isafounder (int id);

/**
 * Find the minimal pedigree connecting the people in the sample list.  The
 * minimal pedigree consists of only those people in the list and their
 * ancestors. Superfluous founders are removed by assigning the child of the
 * superfluous founder pair as a founder.
 *
 * @param nsample
 * @param samplelist
 */
void minimalped(int nsample, int *samplelist);

int findid (int fid);

void pedigree_free (void);

#endif /* PEDIGREE_H */
