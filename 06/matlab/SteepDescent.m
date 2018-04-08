function [inform, x] = SteepDescent(fun, x, sdparams)
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
iter = 0;
alpha_start = 1;

while (iter < sdparams.maxit) && (norm(x.g, Inf) > sdparams.toler * (1 + abs(x.f)))
	d = -x.g;
	[x, alpha, nf, ng] = EBLS(fun, x, d, alpha_start); % x.p, x.g, x.f updated
	iter = iter + 1;
end
inform = struct('status', norm(x.g, Inf) <= sdparams.toler * (1 + abs(x.f)), 'iter', iter);

end
