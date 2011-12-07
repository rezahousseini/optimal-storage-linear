# Storage Optimiziation Problem V2_6
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

/* Storage charge cost */
param pc{i in I};

/* Storage discharge cost */
param pd{i in I};

/* Grid energy import cost */
param pggr{i in I};

/* Fixed energy import cost */
param pgfi{i in I};

/* Free energy import cost */
param pgfr{i in I};

/* Grid energy export cost */
param prgr{i in I};

/* Fixed energy export cost */
param prfi{i in I};

/* Free energy export cost */
param prfr{i in I};

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

/* Maximum free energy delivery change rate */
param Gu{i in I};

/* Maximum free energy request change rate */
param Ru{i in I};

/* Maximum free energy delivery */
param G{i in I};

/* Maximum free energy request */
param R{i in I};

/* Storage charge power */
var uc{i in I}, >= 0, <= min(D, (1/nuc)*gfi[i]/T);

/* Storage discharge power */
var ud{i in I}, >= 0, <= C;

/* Free energy delivery */
var gfr{i in I}, >= 0, <= G[i];

/* Free energy request */
var rfr{i in I}, >= 0, <= R[i];

/* Grid energy delivery */
var ggr{i in I}, >= 0;

/* Grid energy request */
var rgr{i in I}, >= 0;


/* Minimize overall costs */
minimize cost: sum{i in I} (pc[i]*uc[i]+pd[i]*ud[i]+pggr[i]*ggr[i]+prgr[i]*rgr[i]+pgfi[i]*gfi[i]+prfi[i]*rfi[i]+pgfr[i]*gfr[i]+prfr[i]*rfr[i]);

/* Storage degradation constraints */

s.t. Q1: (if Qmax != 0 then Qmax-alphaq*T*(ud[1]-uc[1])/Qmax-betaq*1 else 0), = Qmax;

s.t. Qlimit{k in 2..N}: (if Qmax != 0 then Qmax-alphaq*T*(ud[k]-uc[k])/Qmax-betaq*k else 0), <= (if Qmax != 0 then Qmax-alphaq*T*(ud[k-1]-uc[k-1])/Qmax-betaq*(k-1) else 0);

/* Storage constraints */

/* q_i <= Qmax */
s.t. qub{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), <= (if Qmax != 0 then Qmax-alphaq*T*(ud[i]-uc[i])/Qmax-betaq*i else 0);

/* q_i >= Qmin */
s.t. qlb{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), >= Qmin;

s.t. Dlimit{i in I}: ud[i], <= nud*(nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]))/T;

/* Prosumer node constraints */

s.t. gfr0ub: gfr[1], <= Gu[1];

s.t. rfr0lb: rfr[1], <= Ru[1];

s.t. deltagfrubp{k in 1..N-1}: gfr[k+1], <= gfr[k]+Gu[k+1];

s.t. deltagfrubn{k in 1..N-1}: gfr[k+1], >= gfr[k]-Gu[k+1];

s.t. deltarfrubp{k in 1..N-1}: rfr[k+1], <= rfr[k]+Ru[k+1];

s.t. deltarfrubn{k in 1..N-1}: rfr[k+1], >= rfr[k]-Ru[k+1];

s.t. balance{i in I}: rfi[i]+rfr[i]+rgr[i]-gfi[i]-gfr[i]-ggr[i]+uc[i]-ud[i] = 0;

end;
