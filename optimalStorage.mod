# Storage Optimiziation Problem
#
# Author: Reza Housseini 

/* Set of all time steps */
set I;

/* Set of all storages with Qmax < inf */
set A;

/* Set of all storages with Qmax = inf */
set F;

/* Number of time steps */
param N;

/* Generated energy */
param g{i in I};

/* Requested energy */
param r{i in I};

/* Generated energy cost */
param pg{i in I};

/* Requested energy cost */
param pr{i in I};

/* Storage charge cost */
param pc{a in (A union F), i in I};

/* Storage discharge cost */
param pd{a in (A union F), i in I};

/* Time step */
param T;

/* Maximum storage capacity */
param Qmax{a in (A union F)};

/* Minimum storage capacity */
param Qmin{a in (A union F)};

/* Maximum charge rate */
param C{a in (A union F)};

/* Maximum discharge rate */
param D{a in (A union F)};

/* Efficiency of storage */
param nul{a in (A union F)};

/* Efficiency of charge */
param nuc{a in (A union F)};

/* Efficiency of discharge */
param nud{a in (A union F)};

/* Initial capacity of storage */
param q0{a in (A union F)};

/* Storage charge power */
var uc{a in (A union F), i in I}, >= 0, <= C[a];

/* Storage discharge power */
var ud{a in (A union F), i in I}, >= 0, <= D[a];

/* Real storage level */
var q{a in A, i in I}, >= Qmin[a], <= Qmax[a];

/* Minimize overall costs */
minimize cost: sum{i in I} (pg[i]*g[i]+pr[i]*r[i]+sum{a in (A union F)} (pc[a,i]*uc[a,i]+pd[a,i]*ud[a,i]));

/* Storage constraints */

s.t. qstart{a in A}: q[a,1], = q0[a];

s.t. qnext{a in A, k in 2..N}: q[a,k], = nul[a]*q[a,k-1]+T*(nuc[a]*uc[a,k]-(1/nud[a])*ud[a,k]);

/* Prosumer node constraint */
s.t. balance{i in I}: sum{a in (A union F)}(uc[a,i]-ud[a,i]) = g[i]-r[i];

end;
