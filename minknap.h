// originally part of MINKNAP.C by David Pisinger
// See minknap.c for licensing details

#ifndef MINKNAP_H
#define MINKNAP_H

#define MAXSTATES 400000 

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>


/* ======================================================================
				   macros
   ====================================================================== */

#define SYNC            5      /* when to switch to linear scan in bins */
#define SORTSTACK     200      /* depth of stack used in qsort */
#define MINMED        100      /* find exact median in qsort if larger size */

#define TRUE  1
#define FALSE 0

#define LEFT  1
#define RIGHT 2

#define PARTIATE 1
#define SORTALL  2

#define MAXV (8*sizeof(btype)) /* number of bits in a long integer */
#define PMAX 1                 /* profit of worlds most efficient item  */
#define WMAX 0                 /* weight of worlds most efficient item  */
#define PMIN 0                 /* profit of worlds least efficient item */
#define WMIN 1                 /* weight of worlds least efficient item */

#define DET(a1, a2, b1, b2)    ((a1) * (ptype) (b2) - (a2) * (ptype) (b1))
#define SWAP(a, b)   { register item t; t = *(a); *(a) = *(b); *(b) = t; }
#define DIFF(a,b)              ((int) ((b)-(a)+1))
#define NO(a,p)                ((int) ((p) - (a)->fitem + 1))
#define N(a,p)                 ((int) ((p) - (a)->d.set1))
#define L(x)                   ((long) (x))
#define SZ(a)                  (*(((int *) (a)) - 4) - 1)


/* ======================================================================
				 type declarations
   ====================================================================== */

typedef int           boolean;
typedef long          ntype;   /* number of states/items   */
typedef long          itype;   /* item profits and weights */
typedef long          stype;   /* sum of pofit or weight   */
typedef double        ptype;   /* product type (sufficient precision) */
typedef unsigned long btype;   /* binary representation of solution */

/* item record */
typedef struct irec {
  itype   p;     /* profit */
  itype   w;     /* weight */
  boolean *x;    /* solution variable */
} item;

typedef struct { /* i-stack */
  item  *f;      /* first item in interval */
  item  *l;      /* last item in interval */
} interval;

/* state in dynamic programming */
typedef struct pv {
  stype psum;    /* profit sum */
  stype wsum;    /* weight sum */
  btype vect;    /* solution vector */
} state;

/* set of states */
typedef struct pset {
  ntype size;    /* set size */
  state *fset;   /* first element in set */
  state *lset;   /* last element in set */
  state *set1;   /* first element in array */
  state *setm;   /* last element in array */
} stateset;

typedef struct { /* all problem information */
  ntype    n;               /* number of items         */
  item     *fitem;          /* first item in problem   */
  item     *litem;          /* last item in problem    */
  item     *ftouch;         /* first item considered for reduction */
  item     *ltouch;         /* last item considered for reduction */
  item     *s;              /* current core is [s,t]   */
  item     *t;              /*                         */
  item     *b;              /* break item              */
  item     *fpart;          /* first item returned by partial sort */
  item     *lpart;          /* last item returned by partial sort */
  stype    wfpart;          /* weight sum up to fpart  */
  item     *fsort;          /* first sorted item       */
  item     *lsort;          /* last sorted item        */
  stype    wfsort;          /* weight sum up to fsort  */
  stype    c;               /* current capacity        */
  stype    cstar;           /* origianl capacity       */
  stype    z;               /* current solution        */
  stype    zstar;           /* optimal solution        */
  stype    zwsum;           /* weight sum of zstar     */
  itype    ps, ws, pt, wt;  /* items for deriving bounds */

  btype    vno;             /* current vector number   */
  item *   vitem[MAXV];     /* current last MAXV items */
  item *   ovitem[MAXV];    /* optimal set of items    */
  btype    ovect;           /* optimal solution vector */

  stype    dantzig;         /* dantzig upper bound     */
  stype    ub;              /* global upper bound      */
  stype    psumb;           /* profit sum up to b      */
  stype    wsumb;           /* weight sum up to b      */
  boolean  firsttime;       /* used for restoring x    */
  boolean  welldef;         /* is x welldefined        */
  stateset  d;              /* set of partial vectors  */
  interval *intv1, *intv2;
  interval *intv1b, *intv2b;

  /* debug */
  long     iterates;        /* counters used to obtain specific */
  long     simpreduced;     /* information about the solution process */
  long     pireduced;
  long     pitested;
  long     maxstates;
  long     coresize;
  long     bzcore;
} allinfo;


typedef struct {
	allinfo a;
	item *tab;
	interval *inttab;
} minknap_t;

minknap_t* minknap_prepare(int maxn);
void minknap_free(minknap_t* data);
stype minknap(minknap_t* data, int n, int *p, int *w, int *x, int c);

#endif

