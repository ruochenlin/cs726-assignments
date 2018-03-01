function [x, alpha] = EBLS(fun, x, d, alpha_start)
c1 = 1e-3;
c2 = 0.5;
max_iter = 25;
% x.f = feval(fun, x.p, 1);
% x.g = feval(fun, x.p, 2);
L = 0;
U = Inf;
alpha = alpha_start;
iter = 0;
while iter < max_iter
	step = alpha * d;
	if feval(fun, x.p + step, 1) > x.f + c1 * x.g' * step
		U = alpha;
		alpha = (L + U) / 2;
	elseif feval(fun, x.p + step, 2)' * d < c2 * x.g' * d
		L = alpha;
		if U == Inf
			alpha = 2 * L;
		else
			alpha = (L + U) / 2;
		end
	else
		break;
	end
	iter = iter + 1;
end
x.p = x.p + alpha * d;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
end
