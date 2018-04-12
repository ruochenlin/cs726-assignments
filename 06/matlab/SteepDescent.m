function [inform, x] = SteepDescent(fun, x, sdparams)
global numf; global numg;
numf = 0; numg = 0;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
iter = 0;
alpha_start = 1;

while (iter < sdparams.maxit) && (norm(x.g, Inf) > sdparams.toler * (1 + abs(x.f)))
	% if ~mod(iter,101)
	% 	fprintf('%d: %5.5e, %5.5e\n', iter, x.f, norm(x.g, Inf));
	% end
	d = -x.g;
	[x, alpha, nf, ng] = EBLS(fun, x, d, alpha_start); % x.p, x.g, x.f updated
	iter = iter + 1;
end
inform = struct('status', norm(x.g, Inf) <= sdparams.toler * (1 + abs(x.f)), 'iter', iter);

end
