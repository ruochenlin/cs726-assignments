function [inform, xnew] = LBFGS(fun, x, lbfgsparams)
global numf;
global numg;
numf = 0;
numg = 0;
m_given = lbfgsparams.m;
n = length(x.p);
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
% First step
p = -x.g;
pos_old = x.p;
grad_old = x.g;
[x, steplength] = EBLS(fun, x, p);
s(:, 1) = x.p - pos_old;
y(:, 1) = x.g - grad_old;
rho(:, 1) = 1 / (y(:, 1)' * s(:, 1));
i = 1;
fprintf('%4d: %5.5e, %5.5e\n', i, x.f, norm(x.g));

while norm(x.g) > lbfgsparams.toler * (1 + abs(x.f))
	m = min(m_given, i); % The number of backtrack steps
	index_i = mod(i - 1, m_given) + 1; % The index of info about ith step in arrays y/s/rho 
	gamma = (s(:, index_i)' * y(:, index_i)) / norm(y(:, index_i))^2;
	q = x.g;
	clear alpha;
	for j = i : -1 : i - m + 1
		index_j = mod(j - 1, m_given) + 1; 
		alpha(j - i + m) = rho(:, index_j) * s(:, index_j)' * q; % The subscript of alpha, j-i+m, is the distance between current j and i-m
		q = q - alpha(j - i + m) * y(:, index_j);
	end
	r = gamma * q;
	for j = i - m + 1 : i
		index_j = mod(j - 1, m_given) + 1; 
		beta = rho(:, index_j) * y(:, index_j)' * r;
		r = r + (alpha(j - i + m) - beta) * s(:, index_j);
	end
	% r = H_i * \nabla f_i
	p = -r;
	pos_old = x.p;
	grad_old = x.g;
	[x, steplength] = EBLS(fun, x, p);
	index_i_p_1 = mod(i, m_given) + 1;
	s(:, index_i_p_1) = x.p - pos_old;
	y(:, index_i_p_1) = x.g - grad_old;
	rho(:, index_i_p_1) = 1 / (s(:, index_i_p_1)' * y(:, index_i_p_1));
	i = i + 1;
	fprintf('%4d: %5.5e, %5.5e\n', i, x.f, norm(x.g));
	if i > lbfgsparams.maxit
		break;
	end
end
inform = struct('status', norm(x.g) <= lbfgsparams.toler * (1 + abs(x.f)), 'iter', i);
xnew = x;
end

% No particular use; just to show how to get the correct indeces
function k = get_index(j, m)
	k = mod(j - 1, m) + 1; 
end
