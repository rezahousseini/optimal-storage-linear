function [q, uc, ud, P] = optimalNStorageV1(g, r, p, T, Q, C, D, nul, nuc, nud, q0, cu)
	
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
	% cu: 
	
	N=length(g);
	M=length(Q);
	
	A=zeros((3*N-1)*M,2*N*M);
	B=zeros((3*N-1)*M,1);
	CTYPE=blanks((3*N-1)*M);
	LB=zeros(2*N*M,1);
	UB=zeros(2*N*M,1);
	VARTYPE=blanks(2*N*M);
	COST=zeros(2*N*M,1);
	
	q=zeros(M,N);
	uc=zeros(M,N);
	ud=zeros(M,N);
	
	for m=1:M	
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
				ac(k,l)=nul(m)^(k-l)*nuc(m);
				ad(k,l)=-nul(m)^(k-l)*(1/nud(m));
			endfor
		endfor
	
		% Constraint q_i >= 0
		a(1:N,:)=T*[ac,ad];
		for k=1:N
			b(k)=-nul(m)^k*q0(m);
			ctype(k)="L";
		endfor
		
		% Constraint q_i <= Q
		a(N+1:2*N,:)=T*[ac,ad];
		for k=1:N
			b(N+k)=Q(m)-nul(m)^k*q0(m);
			ctype(N+k)="U";
		endfor
		
		ac=zeros(N);
		ad=zeros(N);
		for k=1:N
			for l=1:k
				ac(k,l)=-nul(m)^(k-l)*nud(m)*nuc(m);
				ad(k,l)=nul(m)^(k-l);
			endfor
			if(k~=N)
				ad(k,k+1)=1;
			endif
		endfor
		
		% Constraint ud_i <= nud*q_i/T
		a(2*N+1:3*N-1,:)=[ac(1:N-1,:),ad(1:N-1,:)];
		for k=1:N-1
			b(2*N+k)=nul(m)^k*nud(m)*q0(m)/T;
			ctype(2*N+k)="U";
		endfor
	
		% Constraint ud_i <= D, uc_i <= min(C, nuc*g_i)
		ubc=zeros(N,1);
		ubd=zeros(N,1);
		for k=1:N
			ubc(k)=min(C(m),nuc(m)*g(k)/T);
			ubd(k)=D(m);
		endfor
		ubd(1)=min(D(m),nud(m)*q0(m)/T);
		lb=zeros(2*N,1);
		ub=[ubc;ubd];
		for k=1:2*N
			vartype(k)="C";
		endfor
		
		% Cost vector
		c=[p,-p]'+cu(m);
		
		A((m-1)*(3*N-1)+1:m*(3*N-1),(m-1)*2*N+1:m*2*N)=a;
		B((m-1)*(3*N-1)+1:m*(3*N-1))=b;
		CTYPE((m-1)*(3*N-1)+1:m*(3*N-1))=ctype;
		LB((m-1)*2*N+1:m*2*N)=lb;
		UB((m-1)*2*N+1:m*2*N)=ub;
		VARTYPE((m-1)*2*N+1:m*2*N)=vartype;
		COST((m-1)*2*N+1:m*2*N)=c;
	endfor
	
	% Solver settings
	s=1;
	param.msglev=1;

	[x, f, status, extra] = glpk (COST, A, B, LB, UB, CTYPE, VARTYPE, s, param); 
	
	% Return values
	for m=1:M
		q(m,1)=q0(m);
		off=(m-1)*2*N;
		uc(m,:)=x(off+1:off+N);
		ud(m,:)=x(off+1+N:off+2*N);
		for k=2:N
			q(m,k)=nul(m)*q(m,k-1)+T*(nuc(m)*uc(m,k-1)-(1/nud(m))*ud(m,k-1));	
		endfor
	endfor
	P=cumsum(p.*r-p.*g+p.*sum(uc,1)-p.*sum(ud,1));
	
endfunction