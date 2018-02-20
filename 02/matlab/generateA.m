function A = generateA(n);
mu=0.01; L=1; kappa=L/mu;
A = randn(n,n); [Q,R]=qr(A);
D=rand(n,1); D=10.^D; Dmin=min(D); Dmax=max(D);
D=(D-Dmin)/(Dmax-Dmin);
D=mu+D*(L-mu);
A=Q'*diag(D)*Q;
end
