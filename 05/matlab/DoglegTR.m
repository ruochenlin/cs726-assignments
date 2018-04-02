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
fprintf(1, ' iter %3d: f=%12.5e, ||Df||=%12.5e, Delta=%7.2e\n', inform.iter, x.f, norm(x.g), delta);
while numf < 100
	inform.iter = inform.iter + 1;
	% Use Dogleg algorithm to solve TR subproblem
	[Q,L] = eig(x.h);
	for i = 1 : length(L)
		if L(i,i) < trparams.delta
			L(i,i) = trparams.delta;
		end
	end
	B = Q * L * Q';
	pu = ((-x.g'*x.g) / (x.g' * B * x.g)) * x.g;
	pb = - B \ x.g;
	norm_pu = norm(pu);
	in_TR = 0;
	if norm_pu > delta
		p = pu * (delta / norm_pu);
	elseif norm(pb) < delta
		p = pb;
		in_TR = 1;
	else
		puTpb = pu' * pb;
		norm_pu_sub_pb = norm(pu - pb);
		alpha = (norm_pu^2 - puTpb + sqrt((puTpb - norm_pu^2)^2 - norm_pu_sub_pb^2 * (norm_pu^2 - delta^2))) / (norm_pu_sub_pb^2);
		p = alpha * pb + (1 - alpha) * pu;
	end
	f_tentative = feval(fun, x.p + p, 1);
	rho = (x.f - f_tentative) / (-0.5 * p' * B * p - x.g' * p);
	if rho < 0.25
		delta = 0.25 * delta;
	else
		if rho > 0.75 && ~in_TR
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
			fprintf(1, ' iter %3d: f=%12.5e, ||Df||=%12.5e, Delta=%7.2e\n', inform.iter, x.f, norm(x.g), delta);
			break;
		end
	end
	fprintf(1, ' iter %3d: f=%12.5e, ||Df||=%12.5e, Delta=%7.2e\n', inform.iter, x.f, norm(x.g), delta);
end
