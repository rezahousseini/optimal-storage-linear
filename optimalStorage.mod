# Storage Optimiziation Problem V2_4
#
# Author: Reza Housseini 

/* Set of all time steps */
set I;

/* Number of all time steps */
param N;

/* */
param gfi{i in I};

/* */
param rfi{i in I};

/* */
param pn{i in I};

/* */
param T;

/* */
param Qmax;

/* */
param Qmin;

/* */
param C;

/* */
param D;

/* */
param nul;

/* */
param nuc;

/* */
param nud;

/* */
param q0;

/* */
param CIst;

/* */
param alphaq;

/* */
param betaq;

/* */
param Gu{i in I};

/* */
param Ru{i in I};

/* */
param G{i in I};

/* */
param R{i in I};

/* Storage charge power */
var uc{i in I}, >= 0, <= min(D, (1/nuc)*gfi[i]/T);

/* Storage discharge power */
var ud{i in I}, >= 0, <= C;

var Q{i in I}, >= Qmin, <= Qmax;

var gfr{i in I}, >= 0, <= G[i];

var rfr{i in I}, >= 0, <= R[i];

/* Minimize overall costs */
minimize cost: sum{i in I} pn[i]*(rfi[i]+rfr[i]-gfi[i]-gfr[i]+T*(uc[i]-ud[i]));

s.t. degrad{i in I}: Q[i], = (if Qmax != 0 then Qmax-alphaq*T*(ud[i]-uc[i])/Qmax-betaq*i else 0);

s.t. Q1: Q[1], = Qmax;

s.t. Qlimit{k in 2..N}: Q[k], <= Q[k-1];

/* q_i <= Qmax */
s.t. qub{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), <= Q[i];

/* q_i >= Qmin */
s.t. qlb{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), >= Qmin;

s.t. Dlimit{i in I}: ud[i], <= nud*(nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]))/T;

s.t. gfr0ub: gfr[1], <= Gu[1];

s.t. rfr0lb: rfr[1], <= Ru[1];

s.t. deltagfrub{k in 1..N-1}: gfr[k+1], <= gfr[k]+Gu[k+1];

s.t. deltarfrub{k in 1..N-1}: rfr[k+1], <= rfr[k]+Ru[k+1];

#s.t. balance{i in I}: rfi[i]+rfr[i]-gfi[i]-gfr[i]+T*(uc[i]-ud[i]) >= 0;

end;
