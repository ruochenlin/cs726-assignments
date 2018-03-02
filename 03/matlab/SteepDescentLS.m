function [inform, x] = SteepDescentLS(fun, x, sdparams)
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
iter = 0;
alpha_start = 1;

while (norm(x.g) > sdparams.toler) && (iter < sdparams.maxit)
	d = -x.g;
	[x, alpha_start] = EBLS(fun, x, d, alpha_start); % x.p, x.g, x.f updated
	if  0. ~= sdparams.eta
		alpha_start = alpha_start * sdparams.eta;
	else
		alpha_start = 1;
	end
	iter = iter + 1;
end
inform = struct('status', norm(x.g) <= sdparams.toler, 'iter', iter);
