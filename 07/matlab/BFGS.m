function [inform, x] = BFGS(fun, x, qnparams)
	global numf;
	global numg;
	numf = 0;
	numg = 0;
	x.f = feval(fun, x.p, 1);
	x.g = feval(fun, x.p, 2);
	m = length(x.p);
	% Use identity as initial approximation of inverse Hessian
	H = eye(m);
	% p is the direction of the next step
	p = -H * x.g;
	pos_old = x.p;
	grad_old = x.g;
	% Perform EBLS along p and update x
	[x, alpha] = EBLS(fun, x, p);
	% Update H_0 with more info
	s = x.p - pos_old;
	y = x.g - grad_old;
	H = (s' * y) / (y' * y) * eye(m);
	iter = 1; fprintf('%4d: %5.5e, %5.5e\n', iter, x.f, norm(x.g));

	while norm(x.g) > qnparams.toler * (1 + abs(x.f))
		% Update H with BFGS algorithm
		rho = 1 / (y' * s);
		HysT = H * y * s';
		H = H + rho * (-HysT - HysT' + (1 + rho*(y'*H*y))*s*s');
		% p is the direction of the next step
		p = -H * x.g;
		pos_old = x.p;
		grad_old = x.g;
		% Perform EBLS along p and update x
		[x, alpha] = EBLS(fun, x, p);
		% Check if maxit has been reached
		iter = iter + 1;
		fprintf('%4d: %5.5e, %5.5e\n', iter, x.f, norm(x.g));

		if iter >= qnparams.maxit
			break;
		end
		s = x.p - pos_old;
		y = x.g - grad_old;
	end
	inform = struct('status', norm(x.g) <= qnparams.toler * (1 + abs(x.f)), 'iter', iter);
end
