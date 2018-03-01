function [inform, x] = SteepDescentBacktrack(fun, x, sdparams)
c1 = 0.001;
beta = 0.5;
alpha = 1;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
iter = 0;
while (norm(x.g) > sdparams.toler) && (iter < sdparams.maxit)
	d = -x.g;
	if 0 ~= sdparams.eta
		alpha = alpha * sdparams.eta;
	else
		alpha = 1;
	end

	while feval(fun, x.p + alpha * d, 1) > (x.f + c1 * x.g' * alpha * d)
		alpha = beta * alpha;
	end
	x.p = x.p + alpha * d;
	x.f = feval(fun, x.p, 1);
	x.g = feval(fun, x.p, 2);
	iter = iter + 1;
end
inform = struct('status', norm(x.g) <= sdparams.toler, 'iter', iter);
