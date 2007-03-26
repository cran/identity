/* pedigree.c - Version 5 */

/* Code for dealing with pedigree information.

   Version 3 differs from version 2 in that the parents within a Person_t
   structure are recorded by their internal id rather than as pointers to
   Person_t structures.

   Version 5 (no version 4) restructures some code to remove the study sample
   code. This information is now included in another file which also includes
   the code to compute the desired quantities.
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
#include "pedigree.h"

#define FOUNDER_ID 0

static Person_t *Pedmem;
static int Npeople;

static void addancestors (int *last_p, int *pedlist, int *inlist);
static int fidcomp (const void *e1, const void *e2);

static idmap_t *Fid2id;
static int *Istruefounder;
static int *Isfounder;
static int *Inminped;
static int *Isinstudy;
static int *Child;
static int *Founderlist;
static int *Minped;
static int *Noffspring;

int
getparent (int whichpar, int indiv)
{
    return Pedmem[indiv].parent[whichpar];
}


void
create_pedigree (int n, int *fid, int *ma, int *pa)
{
    int i, j;
    int found;

    /* the number of people in the pedigree. */
    Npeople = n;

    /* memory allocation. note that we allocate an extra slot since 0 is
     * not used. (founders is 0). */
    Pedmem = Calloc (n + 1, Person_t);
    Fid2id = Calloc (n + 1, idmap_t);
    Istruefounder = Calloc (n + 1, int);
    Isfounder = Calloc (n + 1, int);
    Inminped  = Calloc (n + 1, int);
    Isinstudy = Calloc (n + 1, int);
    Child     = Calloc (n + 1, int);
    Founderlist = Calloc (n + 1, int);
    Minped = Calloc (n + 1, int);
    Noffspring = Calloc (n + 1, int);

    /* zero the first elements of these arrays */
    Pedmem[0].id = 0;
    Pedmem[0].findiv = 0;
    Pedmem[0].parent[0] = Pedmem[0].parent[1] = 0;
    Fid2id[0].id = 0;
    Fid2id[0].fid = 0;
    Istruefounder[0] = 0;
    Isfounder[0] = 0;
    Inminped[0] = 0;
    Isinstudy[0] = 0;
    Child[0] = 0;
    Founderlist[0] = 0;
    Minped[0] = 0;
    Noffspring[0] = 0;

    for (j = 0, i = 1; i <= n; ++j, ++i) {
        Pedmem[i].id = i;
        Pedmem[i].findiv = fid[j];
        Fid2id[i].id = i;
        Fid2id[i].fid = fid[j];

        if (ma[j] == 0 && pa[j] == 0)
            Istruefounder[i] = 1;
        else
            Istruefounder[i] = 0;
    }
    /* sort the mapping hash, so later we can use bisection search */
    qsort (Fid2id, n + 1, sizeof (idmap_t), fidcomp);

    /* Make the parent structure member have the right values. */
    for (i = 0; i < n; ++i) {   /* for each individual */
        if (ma[i] == 0) {
            Pedmem[i+1].parent[0] = 0;
        } else {
            found = 0;
            for (j = 0; j < n; ++j) {
                if (fid[j] == ma[i]) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                error ("Mother of %d not found in pedigree.\n", fid[i]);
            }
            Pedmem[i+1].parent[0] = j + 1;
        }

        if (pa[i] == 0) {
            Pedmem[i+1].parent[1] = 0;
        } else {
            found = 0;
            for (j = 0; j < n; ++j) {
                if (fid[j] == pa[i]) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                error ("Father of %d not found in pedigree.\n", fid[i]);
            }
            Pedmem[i+1].parent[1] = j + 1;
        }
    }
}

int
isafounder (int id)
{
    return (id < 0 || Istruefounder[id] || Isfounder[id]) ? 1 : 0;
}

void
minimalped (int nsample, int *samplelist)
{
    int i, lastidx, nminped, fid, fidx, otherparent;

    Inminped[0] = 1;
    for (i = 1; i <= Npeople; i++) {
        Inminped[i] = 0;
        Isinstudy[i] = 0;
        Child[i] = 0;
        Founderlist[i] = 0;
        Minped[i] = 0;
        Noffspring[i] = 0;
        Isfounder[i] = 0;
    }

    /*
     * Make sure everyone in the sample is in the minimal pedigree
     * and add all of the ancestors for each person.
     */
    for (lastidx = 0, i = 0; i < nsample; i++) {
        Isinstudy[samplelist[i]] = 1;
        if (!Inminped[samplelist[i]]) {
            Inminped[samplelist[i]] = 1;
            Minped[++lastidx] = samplelist[i];
            addancestors(&lastidx, Minped, Inminped);
        }
    }

    /*
     * We want to eliminate founder couples that have only one child
     * and make the child the founder. If the new founder (ie the child)
     * is now part of a founder pair with only one child, repeat the process.
     * To do this, first count the number of offspring everybody has, then
     * go through each founder and check to see if he/she should be removed.
     */
    nminped = lastidx;
    for (fidx = 0, i = 1; i <= nminped; ++i) {
        Noffspring[getparent(0, Minped[i])]++;
        Noffspring[getparent(1, Minped[i])]++;
        Child[getparent(0, Minped[i])] = Minped[i];
        Child[getparent(1, Minped[i])] = Minped[i];
        if (isafounder(Minped[i]))
            Founderlist[fidx++] = Minped[i];
    }
    /*
     * if founder is not in the study sample and has one offspring and
     * the offspring's other parent is a founder not in the study sample
     * with only the offspring as a child, make the offspring a founder
     * and add the offspring to the end of the founder list.
     */
    for (i = 0; (fid = Founderlist[i]) != 0; ++i) {
        if (!Isinstudy[fid] && Noffspring[fid] == 1) {
            otherparent = getparent(0, Child[fid]);
            if (otherparent == fid)
                otherparent = getparent(1, Child[fid]);
            if (isafounder(otherparent) && !Isinstudy[otherparent] &&
                Noffspring[otherparent]==1) {
                Isfounder[Child[fid]] = 1;
                Founderlist[fidx++] = Child[fid];
                Inminped[fid] = 0;
                Inminped[otherparent] = 0;
            }
        }
    }
}

int findid (int fid)
{
   idmap_t target, *result;

   target.fid = fid;
   result = bsearch (&target, Fid2id, Npeople + 1,
                     sizeof (idmap_t), fidcomp);
   if (result) {
       return result->id;
   } else {
       error ("%d not found in pedigree.\n", fid);
       return 0;
   }
}

void
pedigree_free (void)
{
    Free (Pedmem);
    Free (Fid2id);
    Free (Istruefounder);
    Free (Isfounder);
    Free (Inminped);
    Free (Isinstudy);
    Free (Child);
    Free (Founderlist);
    Free (Minped);
    Free (Noffspring);
}

static void
addancestors (int *last_p, int *pedlist, int *inlist)
{
    int idx, mother, father;

    for (idx = *last_p; idx <= *last_p; ++idx) {
        mother = getparent(0, pedlist[idx]);
        if (!inlist[mother]) {
            pedlist[++(*last_p)] = mother;
            inlist[mother] = 1;
        }
        father = getparent(1, pedlist[idx]);
        if (!inlist[father]) {
            pedlist[++(*last_p)] = father;
            inlist[father] = 1;
        }
    }
}

static int
fidcomp (const void *e1, const void *e2)
{
    int v1 = ((idmap_t *) e1)->fid;
    int v2 = ((idmap_t *) e2)->fid;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}
