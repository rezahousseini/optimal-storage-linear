function [q, p, P, uc, ud] = optimal1StorageV2_1(gf, rf, pn, pp, T, Qmax, Qmin, C, D, nul, nuc, nud, q0, CIst, alphaq, betaq)
	
	% INPUT
	% gf: Power generation fix (vector of size N)
	% rf: Power consumption fix (vector of size N)
	% pn: Energy price normal (vector of size N)
	% pp: Energy price power (vector of size N)
	% T: Time period length (scalar)
	% Qmax: Maximum storage capacity (scalar)
	% Qmin: Minimum storage capacity (scalar)
	% C: Maximum charge rate (scalar)
	% D: Maximum discharge rate (scalar)
	% nul: Storage leakage (scalar)
	% nuc: Charge leakage (scalar)
	% nud: Discharge leakage (scalar)
	% q0: Initial storage charge (scalar)
	%
	% OUTPUT
	% q: Storage charge (vector of size N)
	% p: Cost (vector of size N)
	
	N=length(gf);
	
	a=zeros(3*N-1,2*N);
	b=zeros(3*N-1,1);
	ctype=blanks(3*N-1);
	vartype=blanks(2*N);
	lb=zeros(2*N,1);
	ub=zeros(2*N,1);
	c=zeros(2*N,1);
		
	ac=zeros(N);
	ad=zeros(N);
	for k=1:N
		for l=1:k
			ac(k,l)=T*nul^(k-l)*nuc;
			ad(k,l)=-T*nul^(k-l)*(1/nud);
		endfor
	endfor
	
	% Constraint q_i >= Qmin
	a(1:N,:)=[ac,ad];
	for k=1:N
		b(k)=Qmin-nul^k*q0;
		ctype(k)="L";
	endfor
	
	ac=zeros(N);
	ad=zeros(N);
	for k=1:N
		for l=1:k
			ac(k,l)=T*nul^(k-l)*nuc;
			ad(k,l)=-T*nul^(k-l)*(1/nud);
		endfor
		if(k~=N)
			ac(k,k+1)=-alphaq*T/Qmax;
			ad(k,k+1)=alphaq*T/Qmax;
		endif
	endfor
		
	% Constraint q_i <= Q
	a(N+1:2*N,:)=[ac,ad];
	for k=1:N
		b(N+k)=Qmax-nul^k*q0-betaq*(k-1);
		ctype(N+k)="U";
	endfor
		
	ac=zeros(N);
	ad=zeros(N);
	for k=1:N
		for l=1:k
			ac(k,l)=-nul^(k-l)*nud*nuc;
			ad(k,l)=nul^(k-l);
		endfor
		if(k~=N)
			ad(k,k+1)=1;
		endif
	endfor
		
	% Constraint ud_i <= nud*q_i/T
	a(2*N+1:3*N-1,:)=[ac(1:N-1,:),ad(1:N-1,:)];
	for k=1:N-1
		b(2*N+k)=nul^k*nud*q0/T;
		ctype(2*N+k)="U";
	endfor
	
	% Constraint ud_i <= D, uc_i <= min(C, nuc*g_i/T)
	ubc=zeros(N,1);
	ubd=zeros(N,1);
	for k=1:N
		ubc(k)=min(C,nuc*gf(k)/T);
		ubd(k)=D;
	endfor
	ubd(1)=min(D,nud*q0/T);
	lb=zeros(2*N,1);
	ub=[ubc;ubd];
	for k=1:2*N
		vartype(k)="C";
	endfor
	
	% Cost vector
	c=[pn-alphaq*CIst/Qmax^2-pp,-pn+alphaq*CIst/Qmax^2-pp]';
	
	% Solver settings
	s=1;
	param.msglev=1;

	[x, f, status, extra] = glpk (c, a, b, lb, ub, ctype, vartype, s, param);
	
	% Return values
	uc=x(1:N);
	ud=x(N+1:2*N);
	for k=1:N
		p(k)=pn(k)*(rf(k)-gf(k))+(pn(k)-alphaq*CIst/Qmax^2+pp(k))*uc(k)+(-pn(k)+alphaq*CIst/Qmax^2+pp(k))*ud(k)-pp(k)*(uc(k)+ud(k))+pp(k)*(C+D)+betaq*CIst*(k-1)/Qmax;
	endfor
	P=pn*(rf'-gf')+f+sum(betaq*CIst*(0:N)/Qmax)+sum(pp*(C+D));
	q=zeros(1,N+1);
	q(1)=q0;
	for k=2:N+1
		q(k)=nul*q(k-1)+T*(nuc*uc(k-1)-(1/nud)*ud(k-1));
	endfor
	q=q(2:end);
	
endfunction
