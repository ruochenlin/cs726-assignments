function [inform, x] = SteepDescentLS(fun, x, sdparams)
x.g = feval(fun, x.p, 2);
iter = 0;
while (norm(x.g) > sdparams.toler) and (iter < sdparams.maxit)
	x.p = x.p - x.g;
	x.g = feval(fun, x.p, 2);
	iter = iter + 1;
end
inform = struct('status', norm(x.g) <= sdparams.toler, 'iter', iter);
x.f = feval(fun, x.p, 1);
