epsilon = 1e-6;
dim = 100;
globalMinimum = 0;
x_star = zeros(dim, 1);
count_sd = 0;
count_sde = 0;
count_nest = 0;
count_cg = 0;
for i = 1 : 10
	mu=0.01; L=1; kappa=L/mu; n = 100;
	A = randn(n,n); [Q,R]=qr(A);
	D=rand(n,1); D=10.^D; Dmin=min(D); Dmax=max(D);
	D=(D-Dmin)/(Dmax-Dmin);
	D=mu+D*(L-mu);
	A=Q'*diag(D)*Q;
	clear mu; clear L; clear kappa; clear n;
	clear D; clear Dmin;clear Dmax;clear Q;

 	lambda = eig(A);
	L = max(lambda);
	m = min(lambda);
	x0 = randn(dim, 1);
	x_sd(:,:) = x0; x_sde(:,:) = x0; x_nest(:,:) = x0; x_cg(:,:) = x0;
	if 10 == i
		error_sd = []; error_sde = []; error_nest = []; error_cg = [];
		iter_sd = 0; iter_sde = 0; iter_nest = 0; iter_cg = 0;
	end
	% steepest descent
	while (0.5 * x_sd' * A * x_sd - globalMinimum) > epsilon
		if 10 == i
			error_sd(iter_sd + 1, 1) = log10(0.5 * x_sd' * A * x_sd - globalMinimum);
			iter_sd = iter_sd + 1;
		end
		x_sd = x_sd - 1/L * A * x_sd;
		count_sd = count_sd + 1;
	end
	% steepest descent with line search
	while (0.5 * x_sde' * A * x_sde - globalMinimum) > epsilon
		if 10 == i
			error_sde(iter_sde + 1, 1) = log10(0.5 * x_sde' * A * x_sde - globalMinimum);
			iter_sde = iter_sde + 1;
		end
		d = -A * x_sde;
		alpha = -x_sde' * A * d / ( d' * A * d );
		x_sde = x_sde + alpha * d;
		count_sde = count_sde + 1;
	end
	clear d; clear alpha;
	% Nesterov algorithm
	alpha = 1 / L;
	beta = ( sqrt(L) - sqrt(m) ) / ( sqrt(L) + sqrt(m) );
	x_prev(:,:) = x_nest;
	while (0.5 * x_nest' * A * x_nest - globalMinimum) > epsilon
		if 10 == i
			error_nest(iter_nest + 1, 1) = log10(0.5 * x_nest' * A * x_nest - globalMinimum);
			iter_nest = iter_nest + 1;
		end
		x_temp(:,:) = x_nest;
		x_nest = x_nest - alpha * A * ( x_nest + beta * x_nest - x_prev ) + beta * ( x_nest - x_prev );
		x_prev(:,:) = x_temp;
		count_nest = count_nest + 1;
	end
 	clear alpha; clear beta; clear x_temp; clear x_prev;
	% conjugate gradient
	r = A * x_cg;
	p = -r;
	k = 0;
	while (0.5 * x_cg' * A * x_cg - globalMinimum) > epsilon
		if 10 == i
			error_cg(iter_cg + 1, 1) = log10(0.5 * x_cg' * A * x_cg - globalMinimum);
			iter_cg = iter_cg + 1;
		end
		alpha = r' * r / (p' * A * p);
		x_cg = x_cg + alpha * p;
		r_next = r + alpha * A * p;
		beta = r_next' * r_next / ( r' * r );
		p = -r_next + beta * p;
		r = r_next;
		count_cg = count_cg + 1;
	end
	clear r; clear r_next; clear p; clear k; clear alpha; clear beta;
end 

av_sd = count_sd / 10.;
av_sde = count_sde / 10.;
av_nest = count_nest / 10.;
av_cg = count_cg / 10.;

fprintf(1, ' steepest descent - fixed steps : %7.1f\n', av_sd);
fprintf(1, ' steepest descent - exact steps : %7.1f\n', av_sde);
fprintf(1, ' Nesterov                       : %7.1f\n', av_nest);
fprintf(1, ' conjugate gradient             : %7.1f\n', av_cg);
% Generate convergence plot
figure(1);
plot(linspace(0, iter_sd-1, iter_sd), error_sd, linspace(0, iter_sde-1, iter_sde), error_sde, linspace(0, iter_nest-1, iter_nest), error_nest, linspace(0, iter_cg-1, iter_cg), error_cg);
xlabel('k');
ylabel('log10( f(x_k) - f(x*) )');
title('Convergence of four algorithms on f(x) = 1/2*x^TAx');
legend ('SD:const', 'SD:exact', 'Nesterov', 'CG');

x_cg = randn(dim, 1);
r = A * x_cg;
p = -r;
k = 0;
iter_cg_new = 0;
error_cg_new = [];
while iter_cg_new < 120
	error_cg_new(iter_cg_new + 1, 1) = log10(0.5 * x_cg' * A * x_cg - globalMinimum);
	iter_cg_new = iter_cg_new + 1;
	alpha = r' * r / (p' * A * p);
	x_cg = x_cg + alpha * p;
	r_next = r + alpha * A * p;
	beta = r_next' * r_next / ( r' * r );
	p = -r_next + beta * p;
	r = r_next;
end

figure(2);
plot(linspace(0, iter_cg_new - 1, iter_cg_new), error_cg_new);
xlabel('k');
ylabel('log10( f(x_k) - f(x*) )');
title('Convergence of conjugate gradient algorithm on f(x) = 1/2*x^TAx');
clear all;
