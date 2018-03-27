function [inform, x] = DoglegTR(fun, x, trparams)
global numf;
global numg;
global numh;
numf = 0;
numg = 0;
numh = 0;

delta = trparams.Delta0;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
x.h = feval(fun, x.p, 4);
inform = struct('status', 0, 'iter', 0);
for k = 1 : trparams.maxit
	inform.iter = k;
	fprintf(1, ' iter %3d: f=%12.5e, ||Df||=%12.5e, Delta=%7.2e\n', inform.iter, x.f, norm(x.g), delta);
	% Use Dogleg algorithm to solve TR subproblem
	pu = ((-x.g'*x.g) / (x.g' * x.h * x.g)) * x.g;
	pb = - x.h \ x.g;
	norm_pu = norm(pu);
	if norm_pu > delta
		p = pu * (delta * norm_pu);
	elseif norm(pb) <= delta
		p = pb;
	else
		puTpb = pu' * pb;
		norm_pu_sub_pb = norm(pu - pb);
		alpha = (norm_pu^2 - puTpb + sqrt((puTpb - norm_pu^2)^2 - norm_pu_sub_pb^2 * (norm_pu^2 - delta^2))) / (norm_pu_sub_pb^2);
		p = alpha * pb + (1 - alpha) * pu;
	end

	f_tentative = feval(fun, x.p + p, 1);
	rho = (x.f - f_tentative) / (-0.5 * p' * x.h * p - x.g' * p);
	if rho < 0.25
		delta = 0.25 * delta;
	else
		if rho > 0.75 && norm(p) == delta
			delta = min(2 * delta, trparams.hatDelta);
		end
	end
	if rho > trparams.eta
		x.p = x.p + p;
		x.f = f_tentative;
		x.g = feval(fun, x.p, 2);
		x.h = feval(fun, x.p, 4);
		if norm(x.g) <= trparams.toler
			inform.status = 1;
			break;
		end
	end
end
