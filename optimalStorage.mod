# Storage Optimiziation Problem V2_3
#
# Author: Reza Housseini 

/* Set of all time steps */
set I;

/* Number of all time steps */
param N;

/* */
param efi{i in I};

/* */
param pb{i in I};

/* */
param ps{i in I};

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
param Gs{i in I};

/* */
param Rs{i in I};

/* Storage charge power */
var uc{i in I}, >= 0, <= min(D, (1/nuc)*abs(efi[i])/T);

/* Storage discharge power */
var ud{i in I}, >= 0, <= C;

var efr{i in I}, >= Rs[i], <= Gs[i];

var Q{i in I}, >= Qmin, <= Qmax;

/* Minimize overall costs */
minimize cost: sum{i in I} pb[i]*(efi[i]+efr[i]+uc[i]-ud[i]);

s.t. degrad{i in I}: Q[i], = Qmax-alphaq*(ud[i]-uc[i])/(T*Qmax)-betaq*i;

s.t. Q1: Q[1], = Qmax;

s.t. Qlimit{k in 2..N}: Q[k], <= Q[k-1];

/* q_i <= Qmax */
s.t. qub{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), <= Q[i];

/* q_i >= Qmin */
s.t. qlb{i in I}: nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]), >= Qmin;

s.t. Dlimit{i in I}: ud[i], <= nud*(nul^i*q0+T*sum{k in 1..i} nul^(i-1-k)*(nuc*uc[k]-(1/nud)*ud[k]))/T;

#s.t. efr0ub: efr[1], <= Gu[1];

#s.t. efr0lb: efr[1], >= -Ru[1];

s.t. deltaefrub{k in 1..N-1}: efr[k+1], <= efr[k]+Gu[k+1];

s.t. deltaefrlb{k in 1..N-1}: efr[k+1], >= efr[k]-Ru[k+1];

end;
