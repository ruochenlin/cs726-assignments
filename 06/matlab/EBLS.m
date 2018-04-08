function [x, alpha, nf, ng] = EBLS(fun, x, d, alpha_start)
nf = 0; ng = 0;
c1 = 1e-2;
c2 = 0.3;
max_iter = 100;
L = 0;
U = Inf;
alpha = alpha_start;
iter = 0;
while iter < max_iter
	step = alpha * d;
	f_step = feval(fun, x.p + step, 1);
	nf = nf + 1;
	if f_step  > x.f + c1 * x.g' * step
		U = alpha;
		alpha = (L + U) / 2;
	else
		g_step = feval(fun, x.p + step, 2);
		ng = ng + 1;
		if abs(g_step' * d) > abs(c2 * x.g' * d)
			L = alpha;
			if U == Inf
				alpha = 2 * L;
			else
				alpha = (L + U) / 2;
			end
		else
			break;
		end
	end
	iter = iter + 1;
end
x.p = x.p + alpha * d;
x.f = f_step;
x.g = g_step;

end
