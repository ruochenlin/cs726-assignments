function [inform, x] = SteepDescentLS(fun, x, sdparams)
x.g = feval(fun, x.p, 2);
iter = 0;
alpha_start = 1;
while (norm(x.g) > sdparams.toler) && (iter < sdparams.maxit)
	d = -x.g;
	[x, alpha] = EBLS(fun, x, d, alpha_start);
	% d = -x.g;
	x.p = x.p + alpha * d;
	x.g = feval(fun, x.p, 2);
	if  0 ~= sdparams.eta
		alpha_start = alpha * sdparams.eta;
	end
	iter = iter + 1;
end
inform = struct('status', norm(x.g) <= sdparams.toler, 'iter', iter);
x.f = feval(fun, x.p, 1);
