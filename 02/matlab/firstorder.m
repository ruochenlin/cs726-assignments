epsilon = 1e-6;
dim = 100;
globalMinimum = 0;
x_star = zeros(dim, 1);
count_sd = 0;
count_sde = 0;
count_nest = 0;
count_cg = 0;
for i = 1 : 10
	A = generateA(dim); % generate a symmetric positive definite matrix A
	lambda = eig(A);
	L = max(lambda);
	m = min(lambda);
	x0 = randn(dim, 1);
	x_sd(:,:) = x0; x_sde(:,:) = x0; x_nest(:,:) = x0; x_cg(:,:) = x0;
	% steepest descent
	while abs(0.5 * x_sd' * A * x_sd - globalMinimum) > epsilon
		x_sd = x_sd - 1/L * A * x_sd;
		count_sd = count_sd + 1;
	end
	% steepest descent with line search
	while abs(0.5 * x_sde' * A * x_sde - globalMinimum) > epsilon
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
	while abs( 0.5 * x_nest' * A * x_nest - globalMinimum) > epsilon
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
	while ( 0.5 * x_cg' * A * x_cg - globalMinimum) > epsilon
		alpha = r' * r / (p' * A * p);
		x_cg = x_cg + alpha * p;
		r_next = r + alpha * A * p;
		beta = r_next' * r_next / ( r' * r );
		p = -r_next + beta * p;
		r = r_next;
		count_cg = count_cg + 1;
	end
	clear r; clear p; clear k; clear alpha; clear beta;
end 

av_sd = count_sd / 10.;
av_sde = count_sde / 10.;
av_nest = count_nest / 10.;
av_cg = count_cg / 10.;

disp(av_sd);
disp(av_sde);
disp(av_nest);
disp(av_cg);
