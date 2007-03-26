/* nodehash.c */

/*
 * A hash implementation to store and recall kinship probabilities
 * already computed in the nodeprob routine. This version is for
 * computations that are not conditional on genotypes.
 *
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
#include <string.h>
#include <math.h>
#include "ibdgraph.h"
#include "nodehash.h"

typedef  unsigned int  ub4;   /* unsigned 4-byte quantities */
typedef  unsigned char ub1;

static unsigned Ncells;
static unsigned Tablesize;
static unsigned Bsize;
static hash_cell_t *Table;
static unsigned Maxcells;
static unsigned Mask;

static int Nmade;
static int Recycle = 0;
static unsigned Call[15];

static unsigned hashfunc(Ibdgraph_t *node);

static unsigned int
bbhash2(void *key,        /* the key */
        register ub4  length,   /* the length of the key, in ub4s */
        register ub4  initval  /* the previous hash, or an arbitrary value */
        );

#define COLLISION_DEPTH 1
#define areequalnodes(a,b,c,d) ( (memcmp(b,d,sizeof(int)*NGENES)==0) && \
                                 (memcmp(a,c,sizeof(int)*NID)==0) ? 1 : 0 )


void
hashinit (int memsize)
{
    int i;

    for (i=0; i<15; i++)
        Call[i] = 0;

    Nmade = 0;
    Maxcells = ((unsigned) memsize)*1000000/sizeof (hash_cell_t);

    for (Tablesize = 1, i = 1; Tablesize < Maxcells / COLLISION_DEPTH;
         i++, Tablesize<<=1)
        ;
    Bsize = i - 2;
    Tablesize >>= 1;
    /*
     * We want to make sure that the Table size is small enough that there is
     * enough memory available to make a collision list when collisions occur.
     */
    if (Maxcells - Tablesize < Tablesize >> 1) {
        Tablesize >>= 1;
        --Bsize;
    }
    Ncells = Tablesize;
    Mask = ((unsigned)1 << Bsize) - 1;

    Table = Calloc (Tablesize, hash_cell_t);
    if (!Table) {
        error ("Not enough memory for hash table.\n");
    }

    for (i = 0; i < Tablesize; ++i) {
        Table[i].isempty = 1;
        Table[i].next = 0;
        Table[i].prev = 0;
    }
}

#define FNV_INIT ((unsigned)0x811c9dc5)
static unsigned
hashfunc (Ibdgraph_t *node)
{
    int i;
    unsigned hash32;
    static int geneid[NID];
    /*
     * The two component geneid is for the symmetric version
     * where we don't care about the order of the first two ids
     * of the last two ids. If we do care, then just assign each
     * geneid to the values in node->genelist.
     */
    for (i=0; i<NID; i++)
        geneid[i] = node->genelist[i];

    hash32 = bbhash2(geneid, NID, FNV_INIT);
    hash32 = bbhash2(&node->connstate, 1, hash32);

    hash32 &= Mask;

    return hash32;
}

void
hashstore(Ibdgraph_t *node, Probvec_t *prob)
{
    int i;
    hash_cell_t *newcell_p;
    unsigned hv;
    static int geneid[NID];

    /*
     * The two component geneid is for the symmetric version
     * where we don't care about the order of the first two ids
     * of the last two ids. If we do care, then just assign each
     * geneid to the values in node->genelist.
     */
    for (i = 0; i < NID; ++i)
        geneid[i] = node->genelist[i];

    hv = hashfunc (node);
    if (hv > Tablesize){
        error ("hv = %u, too big in hashstore.\n", hv);
    }
    /* Goto entry in hash table and place there, if empty. */
    if (Table[hv].isempty) {
        for (i = 0; i < NCOEF; ++i)
            Table[hv].prob.istate[i] = prob->istate[i];

        Table[hv].cstate = node->connstate;
        for (i = 0; i < NID; ++i)
            Table[hv].nodeid[i] = geneid[i];
        Table[hv].isempty = 0;
        return;
    }
    /*
     * If we haven't hit the max number of cells, create a new one and add
     * it to the beginning of the collision list.
     */
    if (Ncells < Maxcells) {
        newcell_p = (struct hash_cell *) R_alloc (1, sizeof(struct hash_cell));
        if (!newcell_p) {
            error ("Out of memory.\n");
        }
        Ncells++;
        Nmade++;

        newcell_p->next = Table[hv].next;
        for (i=0; i<NCOEF; i++)
            newcell_p->prob.istate[i] = prob->istate[i];
        newcell_p->cstate = node->connstate;
        for (i=0; i<NID; i++)
            newcell_p->nodeid[i] = geneid[i];
        newcell_p->prev = &Table[hv];
        Table[hv].next = newcell_p;
        if (!(newcell_p->next)) /* newcell is the end of collision list */
            Table[hv].prev = newcell_p;
        else
            newcell_p->next->prev = newcell_p;
    } else { /* Are at the max number of cells */
        Recycle = 1;
        /* Place the new cell info at end of list and move it to the front. */
        if (Table[hv].prev && Table[hv].prev != Table[hv].next) {
            /* The case where there are at least two cells in the list. */
            for (i=0; i<NCOEF; i++)
                Table[hv].prev->prob.istate[i] = prob->istate[i];
            /*for (i=0; i<NGENES; i++)
              Table[hv].prev->connectarr[i] = node->connectarr[i];*/
            Table[hv].prev->cstate = node->connstate;
            for (i=0; i<NID; i++)
                Table[hv].prev->nodeid[i] = geneid[i];
            Table[hv].prev->next = Table[hv].next;
            Table[hv].prev->prev->next = 0;
            Table[hv].next = Table[hv].prev;
            Table[hv].prev = Table[hv].prev->prev;
            Table[hv].next->prev = &Table[hv];
            Table[hv].next->next->prev = Table[hv].next;
        } else if (Table[hv].prev) {
            /* Now case where there is only one cell in the collision list. */
            for (i=0; i<NCOEF; i++)
                Table[hv].next->prob.istate[i] = prob->istate[i];
            Table[hv].next->cstate = node->connstate;
            for (i=0; i<NID; i++)
                Table[hv].next->nodeid[i] = geneid[i];
        } else {
            /* No collision list, put everything in the table. */
            for (i=0; i<NCOEF; i++)
                Table[hv].prob.istate[i] = prob->istate[i];
            Table[hv].cstate = node->connstate;
            for (i=0; i<NID; i++)
                Table[hv].nodeid[i] = geneid[i];
        }
    }
}

int hashfind(Ibdgraph_t *node, Probvec_t *prob)
/* Search through hash to find the cell with equivalence classes
   specified in class_p. Returns the probability for the equivalence
   classes. If the search was unsuccessful, returns 0.
*/
{
    int i, connstate=node->connstate;
    int cst = 0;
    int cst2 = 0;
    int *nid, *nid2;
    hash_cell_t *cell_p, *cell_p2;
    unsigned hv;
    static int geneid[NID];
    hash_cell_t *table;
    table = Table;

    /*
     * The two component geneid is for the symmetric version
     * where we don't care about the order of the first two ids
     * of the last two ids. If we do care, then just assign each
     * geneid to the values in node->genelist.
     */
    for (i=0; i<NID; i++)
        geneid[i] = node->genelist[i];

    Call[connstate]++;

    hv = hashfunc(node);
    if (hv > Tablesize) {
        error ("hv = %u too big in hashfind.\n", hv);
    }
    /* Check to see if node is in the table. */
    if (connstate == Table[hv].cstate &&
        (memcmp(geneid,Table[hv].nodeid,sizeof(int)*NID)==0) ) {
        for (i=0; i<NCOEF; i++)
            prob->istate[i] = Table[hv].prob.istate[i];
        return 1;
    }
    else { /* Not in table so check collision list. */
        cell_p = Table[hv].next;
        if (cell_p) {
            cst = cell_p->cstate;
            nid = cell_p->nodeid;
        }
        for ( ; cell_p; cell_p = cell_p2, cst = cst2, nid = nid2) {
            cell_p2 = cell_p->next; //Start getting info for next iteration.
            if (cell_p2) {
                cst2 = cell_p2->cstate;
                nid2 = cell_p2->nodeid;
            }

            if (connstate == cst) {
                if (memcmp(geneid, cell_p->nodeid, sizeof(int)*NID) == 0) {
                    //if (key == nid) {
                    for (i=0; i<NCOEF; i++)
                        prob->istate[i] = cell_p->prob.istate[i];
                    return 1;
                }
            }
        }

    }
    return 0;
}

void
hash_free (void)
{
    Free (Table);
}

/*
  --------------------------------------------------------------------
  mix -- mix 3 32-bit values reversibly.
  For every delta with one or two bit set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
  * If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
  * If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
  mix() was built out of 36 single-cycle latency instructions in a
  structure that could supported 2x parallelism, like so:
  a -= b;
  a -= c; x = (c>>13);
  b -= c; a ^= x;
  b -= a; x = (a<<8);
  c -= a; b ^= x;
  c -= b; x = (b>>13);
  ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
  --------------------------------------------------------------------
*/
#define mix(a,b,c)                              \
    {                                           \
        a -= b; a -= c; a ^= (c>>13);           \
        b -= c; b -= a; b ^= (a<<8);            \
        c -= a; c -= b; c ^= (b>>13);           \
        a -= b; a -= c; a ^= (c>>12);           \
        b -= c; b -= a; b ^= (a<<16);           \
        c -= a; c -= b; c ^= (b>>5);            \
        a -= b; a -= c; a ^= (c>>3);            \
        b -= c; b -= a; b ^= (a<<10);           \
        c -= a; c -= b; c ^= (b>>15);           \
    }

/*
  --------------------------------------------------------------------
  hash() -- hash a variable-length key into a 32-bit value
  k     : the key (the unaligned variable-length array of bytes)
  len   : the length of the key, counting by bytes
  level : can be any 4-byte value
  Returns a 32-bit value.  Every bit of the key affects every bit of
  the return value.  Every 1-bit and 2-bit delta achieves avalanche.
  About 36+6len instructions.

  The best hash table sizes are powers of 2.  There is no need to do
  mod a prime (mod is sooo slow!).  If you need less than 32 bits,
  use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
  In which case, the hash table should have hashsize(10) elements.

  If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

  By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
  code any way you wish, private, educational, or commercial.  It's free.

  See http://burlteburtle.net/bob/hash/evahash.html
  Use for hash table lookup, or anything where one collision in 2^32 is
  acceptable.  Do NOT use for cryptographic purposes.
  --------------------------------------------------------------------
*/

/*
  --------------------------------------------------------------------
  This works on all machines.  hash2() is identical to hash() on
  little-endian machines, except that the length has to be measured
  in ub4s instead of bytes.  It is much faster than hash().  It
  requires
  -- that the key be an array of ub4's, and
  -- that all your machines have the same endianness, and
  -- that the length be the number of ub4's in the key
  --------------------------------------------------------------------
*/
ub4 bbhash2( //k, length, initval)
    void *key,        /* the key */
    register ub4  length,   /* the length of the key, in ub4s */
    register ub4  initval  /* the previous hash, or an arbitrary value */
             )
{
    register ub4 *k = (ub4 *)key;
    register ub4 a,b,c,len;

    /* Set up the internal state */
    len = length;
    a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
    c = initval;           /* the previous hash value */

    /*---------------------------------------- handle most of the key */
    while (len >= 3)
    {
        a += k[0];
        b += k[1];
        c += k[2];
        mix(a,b,c);
        k += 3; len -= 3;
    }

    /*-------------------------------------- handle the last 2 ub4's */
    c += length;
    switch(len)              /* all the case statements fall through */
    {
        /* c is reserved for the length */
    case 2 : b+=k[1];
    case 1 : a+=k[0];
        /* case 0: nothing left to add */
    }
    mix(a,b,c);
    /*-------------------------------------------- report the result */
    return c;
}

/* ======================================================================== */
