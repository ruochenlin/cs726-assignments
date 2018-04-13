function [x, alpha, nf, ng] = EBLS(fun, x, d, alpha_start)
nf = 0; ng = 0;
% EBLS parameters
c1 = 1e-2;
c2 = 0.3;
max_iter = 100;
alpha = alpha_start;

L = 0;
U = Inf;
success = 0;

iter = 0;
while iter < max_iter
	step = alpha * d;
	f_step = feval(fun, x.p + step, 1);
	nf = nf + 1;
	% check the first Wolfe condition
	if f_step  > x.f + c1 * x.g' * step
		% if step is too large
		U = alpha;
		alpha = (L + U) / 2;
	else
		g_step = feval(fun, x.p + step, 2);
		ng = ng + 1;
		% check the second weak Wolfe condition 
		if g_step' * d < c2 * x.g' * d
			% if step is too small
			L = alpha;
			if U == Inf
				alpha = 2 * L;
			else
				alpha = (L + U) / 2;
			end
		else
			% if both Wolfe conditions are satisfied
			success = 1;
			break;
		end
	end
	iter = iter + 1;
end
% update x
x.p = x.p + alpha * d;
if success
	x.f = f_step;
	x.g = g_step;
else
	x.f = feval(fun, x.p, 1);
	nf = nf + 1;
	x.g = feval(fun, x.p, 2);
	ng = ng + 1;
end

end
