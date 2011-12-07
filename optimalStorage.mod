# Storage Optimiziation Problem V1
#
# Author: Reza Housseini 

/* Set of all time steps */
set I;

/* Number of all time steps */
param N;

/* */
param gf{i in I};

/* */
param rf{i in I};

/* */
param pn{i in I};

/* */
param pp{i in I};

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
param Gl{i in I};

/* */
param Gu{i in I};

/* */
param Rl{i in I};

/* */
param Ru{i in I};

/* */
param Gs{i in I};

/* */
param Rs{i in I};

/* Storage charge power */
var uc{i in I}, >= 0, <= min(D, (1/nuc)*gf[i]/T);

/* Storage discharge power */
var ud{i in I}, >= 0, <= C;

var gs{i in I}, >= 0, <= Gs[i];

var rs{i in I}, >= 0, <= Rs[i];

var Q{i in I}, >= Qmin, <= Qmax;

/* Minimize overall costs */
minimize cost: sum{i in I} pn[i]*(rf[i]+rs[i]-(gf[i]+gs[i]-uc[i]+ud[i]));

s.t. degrad{i in I}: Q[i], = Qmax-alphaq*(ud[i]-uc[i])/(T*Qmax)-betaq*i;

s.t. Q1: Q[1], = Qmax;

s.t. Qlimit{k in 2..N}: Q[k], <= Q[k-1];

/* q_i <= Qmax */
s.t. qub{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), <= Q[i];

/* q_i >= Qmin */
s.t. qlb{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), >= Qmin;

s.t. Dlimit{i in I}: ud[i], <= nud*(nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]))/T;

s.t. gs0lb: gs[1], >= -Gl[1];

s.t. gs0ub: gs[1], <= Gu[1];

s.t. rs0lb: rs[1], >= -Rl[1];

s.t. rs0ub: rs[1], <= Ru[1];

s.t. deltagslb{k in 1..N-1}: gs[k+1], >= gs[k]-Gl[k+1];

s.t. deltagsub{k in 1..N-1}: gs[k+1], <= gs[k]+Gu[k+1];

s.t. deltarslb{k in 1..N-1}: rs[k+1], >= rs[k]-Rl[k+1];

s.t. deltarsub{k in 1..N-1}: rs[k+1], <= rs[k]+Ru[k+1];

end;
