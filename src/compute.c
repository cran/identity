/* file: compute.c - Version 1*/

/* All the necessary routines to do the desired computation. In this case,
   we want to compute the identity coefficients for pairs of people.
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

#include <R.h>
#include <Rdefines.h>

#include "ibdgraph.h"
#include "pedigree.h"
#include "nodehash.h"

static void
kin2ident (double *ident, double *istate);

/**
 * Calculate extended identity coefficients.
 *
 * @param fids a vector of integer individual ids.
 * @param ma a vector of integer, mothers' ids.
 * @param pa a vector of integer, fathers' ids.
 * @param pairs a two-column matrix of integers, each line is the pair of ids
 * for which pairwise identify coefficients will be calculated.
 *
 * @return A matrix of 11 columns: id1, id2 and all nine pairwise identity coefs.
 */
SEXP compute_idcoefs (SEXP fids, SEXP ma, SEXP pa, SEXP pairs)
{
   Ibdgraph_t node;
   Probvec_t prob;
   int npairs;                  /* number of pairs */
   int *pair;                   /* pairs[0][i] and pairs[1][i] */
   int pairIndex[2];            /* from user id to internal id */
   SEXP idcoefs;
   double *pcoefs;
   int i, j;

   /* create pedigree structure */
   create_pedigree (LENGTH (fids),
                    INTEGER (fids), INTEGER (ma), INTEGER (pa));

   npairs = INTEGER (GET_DIM (pairs))[1];
   PROTECT (idcoefs = allocMatrix (REALSXP, NIDSTATE, npairs));

   /* initialize hash table, 100 mb */
   hashinit (100);

   /* actual computation */
   pcoefs = REAL (idcoefs);
   pair = INTEGER (pairs);

   for (i = 0; i < npairs; ++i) {
       pairIndex[0] = findid (pair[0]);
       if (pair[1] != pair[0]) {
           pairIndex[1] = findid (pair[1]);
       } else {
           pairIndex[1] = pairIndex[0];
       }
       /* create a minimal pedigree based on two individuals */
       minimalped (2, pairIndex);
       /* main computations */
       ibdgr_init (&node,
                   pairIndex[0], pairIndex[0], pairIndex[1], pairIndex[1]);
       nodeprob (&prob, &node);
       /* save results */
       kin2ident (pcoefs, prob.istate);

       /* moving the pointers */
       pcoefs += NIDSTATE;
       pair += 2;
   }
   /* return */
   UNPROTECT (1);

   pedigree_free ();
   hash_free ();

   return idcoefs;
}

static void
kin2ident (double *ident, double *istate)
{
    int i, j;
    double kincoef[NIDSTATE];
    double kin2idmat[NIDSTATE][NIDSTATE] = {
        {1.00,  0.00, -0.50,  0.00, -0.50,  0.00,  0.50,  0.25,  0.00},
        {0.00,  1.00, -0.50, -1.00, -0.50, -1.00,  0.50,  0.75,  1.00},
        {0.00,  0.00,  2.00,  0.00,  0.00,  0.00, -2.00, -1.00,  0.00},
        {0.00,  0.00,  0.00,  2.00,  0.00,  0.00,  0.00, -1.00, -2.00},
        {0.00,  0.00,  0.00,  0.00,  2.00,  0.00, -2.00, -1.00,  0.00},
        {0.00,  0.00,  0.00,  0.00,  0.00,  2.00,  0.00, -1.00, -2.00},
        {0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  4.00,  0.00,  0.00},
        {0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  4.00,  0.00},
        {0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  4.00} };

    kincoef[0] = istate[0];
    kincoef[1] = istate[1];
    kincoef[2] = istate[2] + istate[3];
    kincoef[3] = istate[4];
    kincoef[4] = istate[5] + istate[6];
    kincoef[5] = istate[7];
    kincoef[6] = istate[8] + istate[9];
    kincoef[7] = istate[10] + istate[11] + istate[12] + istate[13];
    kincoef[8] = istate[14];

    for (i = 0; i < NIDSTATE; ++i) {
        ident[i] = 0;
        for (j = 0; j < NIDSTATE; ++j) {
            ident[i] += kin2idmat[i][j] * kincoef[j];
        }
    }
}
