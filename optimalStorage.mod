# Storage Optimiziation Problem
#
# Author: Reza Housseini 

/* Set of all time steps */
set I;

/* Set of all storages with Qmax < inf */
set F;

/* Set of all storages with Qmax = inf */
set E;

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
param pc{a in (F union E), i in I};

/* Storage discharge cost */
param pd{a in (F union E), i in I};

/* Time step */
param T;

/* Maximum storage capacity */
param Qmax{a in (F union E)};

/* Minimum storage capacity */
param Qmin{a in (F union E)};

/* Maximum charge rate */
param C{a in (F union E)};

/* Maximum discharge rate */
param D{a in (F union E)};

/* Efficiency of storage */
param etal{a in (F union E)};

/* Efficiency of charge */
param etac{a in (F union E)};

/* Efficiency of discharge */
param etad{a in (F union E)};

/* Initial capacity of storage */
param q0{a in (F union E)};

/* Max change of uc */
param DeltaCmax{a in (F union E)};

/* Max change of ud */
param DeltaDmax{a in (F union E)};

/* Storage charge power */
var uc{a in (F union E), i in I}, >= 0, <= C[a];

/* Storage discharge power */
var ud{a in (F union E), i in I}, >= 0, <= D[a];

/* Real storage level */
var q{a in F, i in I}, >= Qmin[a], <= Qmax[a];

/* Minimize overall costs */
minimize cost: sum{i in I} (pg[i]*g[i]+pr[i]*r[i]+sum{a in (F union E)} (pc[a,i]*uc[a,i]+pd[a,i]*ud[a,i]));

/* Storage constraints */

s.t. qstart{a in F}: q[a,1], = q0[a];

s.t. qnext{a in F, k in 2..N}: q[a,k], = etal[a]*q[a,k-1]+T*(etac[a]*uc[a,k]-(1/etad[a])*ud[a,k]);

/*s.t. ucstart{a in (F union E)}: uc[a,1], = 0;*/

s.t. ucnextupper{a in (F union E), k in 2..N}: uc[a,k], <= uc[a,k-1]+DeltaCmax[a];

s.t. ucnextlower{a in (F union E), k in 2..N}: uc[a,k], >= uc[a,k-1]-DeltaCmax[a];

/*s.t. udstart{a in (F union E)}: ud[a,1], = 0;*/

s.t. udnextupper{a in (F union E), k in 2..N}: ud[a,k], <= ud[a,k-1]+DeltaDmax[a];

s.t. udnextlower{a in (F union E), k in 2..N}: ud[a,k], >= ud[a,k-1]-DeltaDmax[a];

/*s.t. uboundupper{a in (F union E), k in 2..N}: (ud[a,k]-uc[a,k])-(ud[a,k-1]-uc[a,k-1]), <= DeltaCmax[a];*/

/*s.t. uboundlower{a in (F union E), k in 2..N}: (ud[a,k]-uc[a,k])-(ud[a,k-1]-uc[a,k-1]), >= -DeltaDmax[a];*/

/* Prosumer node constraint */
/*s.t. balance{i in I}: sum{a in (F union E)}(uc[a,i]-ud[a,i]) = g[i]-r[i];*/
s.t. balanceC{i in I}: sum{a in (F union E)} uc[a,i] = (g[i]-r[i]+abs(g[i]-r[i]))/2;
s.t. balanceD{i in I}: sum{a in (F union E)} ud[a,i] = (r[i]-g[i]+abs(r[i]-g[i]))/2;

end;
