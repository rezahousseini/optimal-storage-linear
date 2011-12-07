function [q, P] = optimal1StorageV1(g, r, p, T, Q, C, D, nul, nuc, nud, q0)
	
	% g:
	% p:
	% T:
	% Q:
	% C:
	% D:
	% nul:
	% nuc:
	% nud:
	% q0:
	% [q_N: 
	% P: Profit
	% q: Storage charge
	
	N=length(g);
	
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
			ac(k,l)=nul^(k-l)*nuc;
			ad(k,l)=-nul^(k-l)*(1/nud);
		endfor
	endfor
	
	% Constraint q_i >= 0
	a(1:N,:)=T*[ac,ad];
	for k=1:N
		b(k)=-nul^k*q0;
		ctype(k)="L";
	endfor
		
	% Constraint q_i <= Q
	a(N+1:2*N,:)=T*[ac,ad];
	for k=1:N
		b(N+k)=Q-nul^k*q0;
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
	
	% Constraint ud_i <= D, uc_i <= min(C, (1/nuc)*g_i/T)
	ubc=zeros(N,1);
	ubd=zeros(N,1);
	for k=1:N
		ubc(k)=min(C,(1/nuc)*g(k)/T);
		ubd(k)=D;
	endfor
	ubd(1)=min(D,nud*q0/T);
	lb=zeros(2*N,1);
	ub=[ubc;ubd];
	for k=1:2*N
		vartype(k)="C";
	endfor
		
	% Cost vector
	c=[p,-p]';
	
	% Solver settings
	s=1;
	param.msglev=1;

	[x, f, status, extra] = glpk (c, a, b, lb, ub, ctype, vartype, s, param); 
	
	% Return values
	P=p*(r'-g')+f;
	uc=x(1:N);
	ud=x(N+1:2*N);
	q=zeros(1,N+1);
	q(1)=q0;
	for k=2:N+1
		q(k)=nul*q(k-1)+T*(nuc*uc(k-1)-(1/nud)*ud(k-1));	
	endfor
	q=q(2:end);
	
endfunction
