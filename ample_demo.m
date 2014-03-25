% ample_demo.m
%
% A demonstration test-case for ample.m. Will
% calculate the factorization of the posterior 
% distribution of an iid Gauss-Bernoulli
% signal sampled with an iid Gaussian random matrix 
% given the set of observations.

%% Demo Parameters
N = 2^12;				% Signal dimensionality
subrate  = 0.50;		% Ratio of M/N (percent of dim. reduction)
sparsity = 0.25;		% Percent of signal which is non-zero
gb_mean = 0.5;			% Mean of GB signal prior
gb_var  = 1;			% Variance of GB signal prior
delta   = 1e-8;			% iid AWGN variance 
M = round(N*subrate);	% Number of measurements
K = round(sparsity*N);	% Number of non-zeros

%% Generate Problem
A = randn(M,N) ./ sqrt(N);		  % A random iid projector		
x  = sqrt(gb_var).*randn(N,1);	  % Generate gaussian part of signal...
rp = randperm(N);  				  % Get random nonzero locations...
z = rp(K+1:end); 
nz = rp(1:K);
x(z) = 0;						  % Set the zeros to make the signal sparse.
y = A*x + sqrt(delta)*randn(M,1); % Calculate noisy measurements

%% Solve with ample-GB
fprintf('Running ample-GB...\n');
[a_gb,c_gb,history_gb] = ample(	A,y,@prior_gb,...
					     		[gb_mean, gb_var, sparsity],...
					  			x);


%% Solve with ample-L1
fprintf('Running ample-L1...\n');
[a_l1,c_l1,history_l1] = ample(	A,y,@prior_l1sparse,...
					     		[-5 5],...
					  			x);

%% Reporting
fprintf('-----------------------\n');
mse_gb = norm(a_gb - x).^2./N;
fprintf('(GB) Final MSE : %0.2e\n',mse_gb);
mse_l1 = norm(a_l1 - x).^2./N;
fprintf('(L1) Final MSE : %0.2e\n',mse_l1);

%% Display
loc = 1:1:N;
figure(1); clf;
% Plot the ample-GB results
subplot(2,1,1);
	hold on;
		stem(loc(nz),x(nz),'-bo','DisplayName','Original');
		stem(loc,a_gb,'-rx','MarkerSize',3,'DisplayName','Recovered Means');
	hold off;
	grid on; box on;
	xlabel(sprintf('MSE = %0.2e',mse_gb));
	title('Recovery with true GB prior');
	legend('Location','EastOutside');
	axis([1 N -3 3]);
% Plot the ample-L1 results
subplot(2,1,2);
	hold on;
		stem(loc(nz),x(nz),'-bo','DisplayName','Original');
		stem(loc,a_l1,'-rx','MarkerSize',3,'DisplayName','Recovered Means');
	hold off;
	grid on; box on;
	xlabel(sprintf('MSE = %0.2e',mse_l1));
	title('Recovery with L1-Sparse prior');
	legend('Location','EastOutside');
	axis([1 N -3 3]);