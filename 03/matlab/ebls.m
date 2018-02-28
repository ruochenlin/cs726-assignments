function [x, alpha] = EBLS(fun, x, d, alpha_start)
c1 = 1e-3;
c2 = 0.5;
max_iter = 25;
L = 0;
U = Inf;
alpha = alpha_start;
iter = 1;
while iter <= max_iter
	x.f = feval(fun, x.p, 1);
	x.g = feval(fun, x.p, 2);
	step = -alpha * x.g;
	if feval(fun, x.p + step, 1) > x.f + c1 * x.g' * step
		U = alpha;
		alpha = (L + U) / 2;
	else 
		if feval(fun, x.p + step, 2)'*step < c2 * x.g' * step
			L = alpha;
			if U = Inf
				alpha = 2 * L;
			else
				alpha = (L + U) / 2;
			end
		else
			break;
		end
	end
end
