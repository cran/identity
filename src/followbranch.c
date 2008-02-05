/* followbranch.c - Version 7 */

/* Recursive routine for replacing an individual with his/her
   parents within a list of individuals. This action must never
   be taken on a person who has descendants in the list. That is,
   the recursion must work from the bottom of the pedigree (youngest
   individuals) up. The recursion stops when only founders are left
   in the list of people. When a person appears more than once in the
   list, all instances of that person must be replaced by one of the
   parents simultaneously. This results in 2^m possible branches to
   follow where m is the number of times the indivdual appears in the
   list. In cases where multiple instances of the same person get
   replaced by the same parent, then all those instances of the parent
   become "connected." All connected instances of a person must
   get replaced by the same parent. The actual number of branches
   to follow in this case becomes 2^k, where k is the number of
   disconnected groups.
*/
/*
  Copyright (C) 2004, 2005 Mark Abney

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 * As of May 13, 2004 followbranch1.c is verified to work properly.
 * Subsequent versions will focus on greater speed by specializing
 * to the case of four alleles with symmetry within the first two positions
 * and within the last two positions, memory allocation optimizations
 * within the recursive algorithm, minimizing unnecessary data copying/transfer
 * between functions, implementing a hash table to reduce recomputation for
 * other pairs, and possibly other forms of speed up.
 */

/*
 * Version 2 will implemention various optimization by reducing the amount
 * of memory allocations and data copying.
 *
 * Version 3 includes storing computed probabilities in a hash table,
 * and recalling these probabilities when needed. Appear to work properly
 * as of May 24, 2004. Requires pedigree.c, ibdgraph.c, nodehash.c, and
 * tested with maintest2.c and identkin.c.
 *
 * Version 3.1 uses a symmetric version of nodehash such that the order
 * of individuals in the first two entries of genenode->genelist is not
 * relevant, and the order of people in the last two entries of genelist
 * are also not important.
 *
 * Version 4 is developed from 3.1
 *
 * Version 5 is essentially identical to version 4.1.
 *
 * Version 7 reduces the number of branches to follow, or look up in the
 * hash table by using the symmetry in nodes such as (6 1 6 2). In this case
 * the only branches followed are (3 1 3 2), (3 1 4 2), (4 1 4 2). The
 * probabilities for node (4 1 3 2) are computed as a permutation of the
 * second node. This is also done for nodes such as (6 6 6 1) and
 * (6 6 6 6) as long as none of the 6's are connected.
 *
 * Version 8 reimplements the symmetry calculations that were removed for
 * version 7. Now, the genelist in the node is reordered from largest to
 * smallest id, with the connection array also reordered. The resulting
 * permutation is inverted prior to returning the probability vector so that
 * the proper values are passed back to the calling function.
 */

#include <R.h>
#include <Rdefines.h>

#include "pedigree.h"
#include "ibdgraph.h"
#include "nodehash.h"

/* #define DEBUG 1 */

/**
 * Given the probability vector for the connection states of four genes,
 * return the probability vector associated with the states that result when
 * genes idx1 and idx2 are switched.
 *
 * @param idx1
 * @param idx2
 * @param prob
 * @param probnew
 *
 */
static void
pswap (int idx1, int idx2, const Probvec_t *prob, Probvec_t *probnew);

/**
 * Whether to do the symmetrized computation (e.g. only follow three branches
 * when the node is (6 1 6 2)).
 *
 * @param node
 * @param slot
 * @param ncop
 *
 * @return
 */
static int
dosymcomp (Ibdgraph_t *node, int *slot, int ncop);

/**
 * Compute the probabilities of the branches not taken in the symmetric
 * case and add them to the probability of the node that was computed.
 *
 * @param branch
 * @param ncop
 * @param slot
 * @param prob
 */
static void
symprobs (int branch, int ncop, int *slot, Probvec_t *prob);

/**
 * Sort the genelist based on (1) id number, (2) number of connections,
 * (3) element number in the connection array.
 *
 * @param node
 * @param perm
 */
static void
sortnode (Ibdgraph_t *node, int *perm);

static int
lessthan(int a1, int a2, int a3, int v1, int v2, int v3);

/**
 * Permute the probability vector so it corresponds to the probability
 * of the node before it was sorted.
 *
 * @param prob
 * @param perm
 */
static void
unsortprob(Probvec_t *prob, int *perm);

static int
arerepeatedfounders (Ibdgraph_t *node, int *slot);

static int
onlyfounders(int n, int *ids);

static void
initpvec (Probvec_t *pvec);

static void
add2pvec (Probvec_t *orig, Probvec_t *x);

static void
mult2pvec (double coef, Probvec_t *x);

static void
newboundcond(Ibdgraph_t *node, Probvec_t *prob);

static void
connected2parent(Ibdgraph_t *node, int prnt, int idx);

static void
founder2parent(Ibdgraph_t *node, int prnt, int idx);

static void
connect_id (Ibdgraph_t *node, int idx1, int idx2);

static int
connectstate(int *connect);

#define areconnected(c,i,j) ( (c)[i] == (c)[j] )
//int areconnected(int *connectar, int idx1, int indx2);
#define SWAP(a,b,tmp) {(tmp)=(a); (a)=(b); (b)=(tmp);}

void ibdgr_init (Ibdgraph_t *pnode,
                 int id1, int id2, int id3, int id4)
{
    int i;

    pnode->genelist[0] = id1;
    pnode->genelist[1] = id2;
    pnode->genelist[2] = id3;
    pnode->genelist[3] = id4;

    for (i = 0; i < NGENES; ++i) {
        pnode->connectarr[i] = i;
    }
}

void
nodeprob (Probvec_t *pprob, Ibdgraph_t *genenode)
{
    int i, j, k, number, unconnected;
    double coef;

    int parent[NGENES];         /* A binary number with each element in the
                                   array being a 0 or 1. This number specifies
                                   which parent to branch to for each
                                   unconnected copy of the branching person. */

    int index[NGENES];          /* index[0] is the index of the largest
                                   id in geneset.genelist. index[1] is the
                                   second largest, etc. */
    int slot[NGENES];          /* slot[0] is the index of the person to branch
                                  on in the genelist. slot[1] is the index of
                                  the second instance of that person in the list
                                  that is not connected to the first. */
    int ncop; /* number of unconnected copies of branching person. */
    int nrep = 0; /* number of times a founder is repeated (=ncop-1) */
    Ibdgraph_t newnode;         /* The new (post branching) list is placed
                                   here and passed on to do another branching.*/
    Probvec_t tmppvec;

    int symbranch;              /* Whether to do symmetric branching. */
    int nbranch;                /* number of branches to take. */

    initpvec (pprob);

    /*
     * We need to sort the node here so that any subsequent calls to hashfind
     * will find the node (hashfind assumes the node is already sorted). The
     * probability vector will need to be permuted to match the unsorted node
     * after storing in the hash table, but before returning the probability.
     */
    sortnode (genenode, index);
    genenode->connstate = connectstate (genenode->connectarr);

    /*
     * The boundary conditions apply when all individuals in the node
     * are founders.
     */
    if (onlyfounders (NGENES, genenode->genelist)) {
        nrep = arerepeatedfounders (genenode, slot);
        if (!nrep) {
            /* No more branching is possible, so compute probability of
             * the configuration defined in the connection array of genenode.*/
            newboundcond(genenode, pprob);
            // We're not bothering to store the boundary condition probabilities.
            unsortprob(pprob, index);
            return;
        }
        /* There are only founders in the node, but some are repeated. */
        ncop = nrep + 1;

        /*
         * First check if the probability for this node has been computed,
         * and only if it has not do the computation.
         */
        if (hashfind(genenode, pprob)) {
            unsortprob(pprob, index);
            return;
        }
    } else { // There are some nonfounders in the node.
        /*
         * First check if the probability for this node has been computed,
         * and only if it has not do the computation.
         */
        if (hashfind(genenode, pprob)) {
            unsortprob(pprob, index);
            return;
        }

        /*
         * The algorithm proceeds by first sorting the indices of the
         * individuals in the genelist, to find the person to branch
         * on. Then for each unconnected copy of this person in the
         * genelist, replace the person with one of the parents. The branch
         * that we will follow is represented by a binary number where each
         * digit in the binary number represents a replacement of the person
         * with either the mother (0) or father(1). We insure following all
         * branches by looping through all possible binary numbers.
         */
        for (i = 0; isafounder (genenode->genelist[i]) && i<NGENES; i++)
            ; //  for loop is deliberately empty
        slot[0] = i; // slot[0] has index of largest id that is not a founder
        ncop = 1;
        ++i;
        for (; i < NGENES && genenode->genelist[i] == genenode->genelist[i-1];
             ++i) {
            if (!areconnected(genenode->connectarr, i-1, i)) {
                for (unconnected = 1,j = i-1; unconnected && j >= 0; --j) {
                    if (areconnected(genenode->connectarr, j, i))
                        unconnected = 0;
                }
                if (unconnected) {
                    slot[ncop] = i;
                    ncop++;
                }
            }
        }
    }

    /*
     * Now we want to loop through all the branches that need to be
     * followed. Making sure to connect unconnected genes that branch up
     * to the same parent, and replacing all connected genes by the same
     * parent.
     */
    coef = 1.0 / (1 << ncop);
    if (dosymcomp (genenode, slot, ncop)) {
        symbranch = 1;
        nbranch = ncop+ 1 ;
    } else {
        symbranch = 0;
        nbranch = 1 << ncop;
    }
    for (i = 0; i < nbranch; i++) { //Loop over branches
        newnode = *genenode;

        /* Succesively peel off a 0 or 1 for each copy of the gene in the list.*/
        if (symbranch) {
            number = (1 <<i ) - 1;
        } else {
            number = i;
        }

        for (j = 0; j < ncop; ++j) { //Replace each copy of person with a parent
            parent[j] = number & 1;
            number >>= 1;
            if (nrep) // nrep nonzero only when branching on a repeated founder
                founder2parent (&newnode, parent[j], slot[j]);
            else
                connected2parent (&newnode, parent[j], slot[j]);
            if (j > 0) {
                /* Move k to next gene copy going to same parent, connect.*/
                for (k = j - 1; k>=0 && parent[k] != parent[j]; k--) ;
                if (k >= 0)
                    connect_id (&newnode, slot[k], slot[j]);
            }
        }
        nodeprob (&tmppvec, &newnode);
        if (symbranch && i != 0 && i != nbranch-1) {
            symprobs (i, ncop, slot, &tmppvec);
        }
        add2pvec (pprob, &tmppvec);
    }

    mult2pvec (coef, pprob);
    hashstore (genenode, pprob);
    unsortprob (pprob, index);
}

static int
arerepeatedfounders (Ibdgraph_t *node, int *slot)
{
    int found = 0;
    int nrep = 0;
    int i, j, k, unconnected, id1;

    for (i = 0; i < NGENES && !found; ++i) {
        id1 = node->genelist[i];
        slot[0] = i;
        for (j = i + 1; j < NGENES; ++j) {
            if (node->genelist[j] == id1 &&
                !areconnected(node->connectarr,i,j)) {
                for (unconnected = 1, k = j-1; unconnected && k >= 0; --k) {
                    if (areconnected(node->connectarr, j, k))
                        unconnected = 0;
                }
                if (unconnected) {
                    found = 1;
                    nrep++;
                    slot[nrep] = j;
                }
            }
        }
    }
    return nrep;
}

static int
onlyfounders (int n, int *ids)
{
    int nfound = 0;
    int all = 1;
    int i;

    for (i = 0; i < n && all; ++i) {
        if (!isafounder(ids[i])) {
            all = 0;
        }
    }
    return all;
}

static int
dosymcomp (Ibdgraph_t *node, int *slot, int ncop)
{
    int i, j, idx;

    for (i = 0; i < ncop; ++i) {
        idx = slot[i];
        for (j = idx + 1; j < NGENES; ++j)
            if (areconnected(node->connectarr, idx, j))
                return 0;
    }
    return 1;
}

static void
symprobs (int branch, int ncop, int *slot, Probvec_t *prob)
{
    static Probvec_t pr1, pr2, tempprob;

    switch (ncop) {
    case 2:
        pswap (slot[0], slot[1], prob, &tempprob);
        break;
    case 3:
        switch (branch) {
        case 1: // origprob = prob(branch 100)
            pswap (slot[0], slot[1], prob, &tempprob); // branch 010
            break;
        case 2: // origprob = prob(branch 110)
            pswap (slot[1], slot[2], prob, &tempprob); // branch 101
            break;
        }
        pswap (slot[0], slot[2], prob, &pr1); // branch 001/011
        add2pvec (&tempprob, &pr1);
        break;
    case 4:
        switch (branch) {
        case 1: // origprob = prob(branch 1000)
            pswap (slot[0], slot[1], prob, &tempprob); // branch 0100
            pswap (slot[0], slot[2], prob, &pr1); // branch 0010
            add2pvec (&tempprob, &pr1);
            pswap (slot[0], slot[3], prob, &pr1); // branch 0001
            add2pvec (&tempprob, &pr1);
            break;
        case 2: // origprob = prob(branch 1100)
            pswap (slot[1], slot[2], prob, &tempprob); // branch 1010
            pswap (slot[1], slot[3], prob, &pr1); // branch 1001
            add2pvec (&tempprob, &pr1);

            pswap (slot[0], slot[2], prob, &pr1); // branch 0110
            add2pvec (&tempprob, &pr1);

            pswap (slot[0], slot[3], prob, &pr1); // branch 0101
            add2pvec (&tempprob, &pr1);
            pswap (slot[1], slot[2], &pr1, &pr2); // branch 0011
            add2pvec (&tempprob, &pr2);
            break;
        case 3: // origprob = prob(branch 1110)
            pswap (slot[2], slot[3], prob, &tempprob); // branch 1101

            pswap (slot[1], slot[3], prob, &pr1); // branch 1011
            add2pvec (&tempprob, &pr1);
            pswap (slot[0], slot[3], prob, &pr1); // branch 0111
            add2pvec (&tempprob, &pr1);
            break;
        }
        break;
    }
    add2pvec (prob, &tempprob);
}

static void
sortnode (Ibdgraph_t *node, int *perm)
{
    int *list = node->genelist;
    int *carr = node->connectarr;
    int count[4] = {0, 0, 0, 0};
    int i, j;
    int val1, val2, val3;

    // First need to count the number of connections each element has.
    for (i = 0; i < NGENES; ++i)
        count[carr[i]]++;

    perm[0] = 0;
    for (i=1; i<NGENES; i++) {
        val1 = list[i];
        val2 = count[carr[i]];
        val3 = carr[i];
        for (j = i;
             j > 0 && lessthan (list[j-1], count[carr[j-1]], carr[j-1],
                                val1, val2, val3);
             j--) {
            perm[j] = perm[j-1];
            list[j] = list[j-1];
            carr[j] = carr[j-1];
        }
        perm[j] = i;
        list[j] = val1;
        carr[j] = val3;
    }
}

static int
lessthan (int a1, int a2, int a3, int v1, int v2, int v3)
{
    if (a1 < v1) {
        return 1;
    } else if (a1 > v1) {
        return 0;
    } else if (a2 < v2) {
        return 1;
    } else if (a2 > v2) {
        return 0;
    } else if (a3 > v3) {
        return 1; // Greater than to stop extra swapping.
    } else {
        return 0;
    }
}

static void
unsortprob (Probvec_t *prob, int *perm)
{
    int i, pidx, tmp;
    int done;
    do {
        done = 1;
        for (i = 0; i < NGENES; ++i) {
            pidx = perm[i];
            if (pidx != i) {
                done = 0;
                pswap (i, pidx, prob, prob);
                SWAP(perm[i], perm[pidx],tmp);
            }
        }
    } while (!done);
}

static void
pswap (int idx1, int idx2, const Probvec_t *prob, Probvec_t *newprob)
{
    const double *pvec_old = prob->istate;
    double *pvec = newprob->istate;
    double tmp;
    int itmp;
    int i;

    if (idx1 > idx2) SWAP(idx1, idx2, itmp);

    for (i = 0; i < NCOEF; ++i) {
        pvec[i] = pvec_old[i];
    }

    switch (idx1) {
    case 0:
        switch (idx2) {
        case 1:
            SWAP(pvec[5], pvec[6],tmp);
            SWAP(pvec[8], pvec[9],tmp);
            SWAP(pvec[10], pvec[12],tmp);
            SWAP(pvec[11], pvec[13],tmp);
            break;
        case 2:
            SWAP(pvec[1], pvec[9],tmp);
            SWAP(pvec[3], pvec[6],tmp);
            SWAP(pvec[4], pvec[12],tmp);
            SWAP(pvec[7], pvec[11],tmp);
            break;
        case 3:
            SWAP(pvec[1], pvec[8],tmp);
            SWAP(pvec[2], pvec[6],tmp);
            SWAP(pvec[4], pvec[13],tmp);
            SWAP(pvec[7], pvec[10],tmp);
            break;
        }
        break;
    case 1:
        switch (idx2) {
        case 2:
            SWAP(pvec[1], pvec[8],tmp);
            SWAP(pvec[3], pvec[5],tmp);
            SWAP(pvec[4], pvec[10],tmp);
            SWAP(pvec[7], pvec[13],tmp);
            break;
        case 3:
            SWAP(pvec[1], pvec[9],tmp);
            SWAP(pvec[2], pvec[5],tmp);
            SWAP(pvec[4], pvec[11],tmp);
            SWAP(pvec[7], pvec[12],tmp);
            break;
        }
        break;
    case 2:
        SWAP(pvec[2], pvec[3],tmp);
        SWAP(pvec[8], pvec[9],tmp);
        SWAP(pvec[10], pvec[11],tmp);
        SWAP(pvec[12], pvec[13],tmp);
        break;
    }
}


static void
initpvec (Probvec_t *pvec)
{
    int i;
    for (i = 0; i < NCOEF; ++i) {
        pvec->istate[i] = 0.0;
    }
}

static void
add2pvec(Probvec_t *orig, Probvec_t *x)
{
    int i;
    for (i = 0; i < NCOEF; ++i) {
        orig->istate[i] += x->istate[i];
    }
}

static void
mult2pvec (double coef, Probvec_t *x)
{
    int i;
    for (i = 0; i < NCOEF; ++i) {
        x->istate[i] *= coef;
    }
}


static void
newboundcond (Ibdgraph_t *node, Probvec_t *prob)
{
   int state = connectstate(node->connectarr);
   prob->istate[state] = 1;
}

static void
connected2parent (Ibdgraph_t *node, int prnt, int idx)
{
    int i, id;

    id = node->genelist[idx];
    for (i = 0; i < NGENES; i++) {
        if (areconnected(node->connectarr, idx, i) && node->genelist[i] == id)
            node->genelist[i] = getparent(prnt, node->genelist[i]);
    }
}

static void
founder2parent (Ibdgraph_t *node, int prnt, int idx)
{
    int i, id;

    id = node->genelist[idx];
    for (i = 0; i < NGENES; ++i) {
        if (areconnected(node->connectarr, idx, i) && node->genelist[i] == id)
            if (prnt)
                node->genelist[i] = - node->genelist[i];
    }
}

static void
connect_id (Ibdgraph_t *node, int idx1, int idx2)
{
    int *conn_state = node->connectarr;
    int i, state;

    if (conn_state[idx1] != conn_state[idx2]) {
        for (state = conn_state[idx2], i = 0; i < NGENES; ++i) {
            if (conn_state[i] == state) {
                conn_state[i] = conn_state[idx1];
            }
        }
    }
}

static int
connectstate(int *connect)
{
    static int state[64] = { 14,  7, 13, -1, 12, -1, -1,  6, //0-7
                             11, -1, -1, -1,  9, -1, -1, -1, //8-15
                             10, -1,  8, -1, -1, -1, -1, -1, //16-23
                             -1,  5, -1, -1, -1, -1, -1, -1, //24-31
                             4,  1, -1, -1, -1, -1, -1, -1, //32-39
                             -1, -1,  3, -1, -1, -1, -1, -1, //40-47
                             -1, -1, -1, -1,  2, -1, -1, -1, //48-55
                             -1, -1, -1, -1, -1, -1, -1,  0 }; //56-63

    int ab = connect[0] == connect[1];
    int ac = connect[0] == connect[2];
    int ad = connect[0] == connect[3];
    int bc = connect[1] == connect[2];
    int bd = connect[1] == connect[3];
    int cd = connect[2] == connect[3];

    int idx = (ab<<5) + (ac<<4) + (ad<<3) + (bc<<2) + (bd<<1) + cd;

    return state[idx];
}
