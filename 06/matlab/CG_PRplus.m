function [inform, x] = CG_PRplus(fun, x, nonCGparams)
global numf; global numg;
numf = 0; numg = 0;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
alpha_start = 1;
p = -x.g;
iter = 0;
while (iter < nonCGparams.maxit) && (norm(x.g, Inf) > nonCGparams.toler * (1 + abs(x.f)))
	% if ~mod(iter,10)
	% 	fprintf('%d: %5.5e, %5.5e\n', iter, x.f, norm(x.g, Inf));
	% end
	
	g_old = x.g;
	[x, alpha, nf, ng] = EBLS(fun, x, p, alpha_start);
	beta = x.g' * (x.g - g_old) / (g_old' * g_old);
	if beta > 0
		p = -x.g + beta * p;
	else
		p = -x.g;
	end
	
	iter = iter + 1;
end
inform = struct('status', norm(x.g, Inf) <= nonCGparams.toler * (1 + abs(x.f)), 'iter', iter);

end 
