% ample_gb_demo.m
%
% A demonstration test-case for ample.m. Will
% calculate the factorization of the posterior 
% distribution of an iid Gauss-Bernoulli
% signal sampled with an iid Gaussian random matrix 
% given the set of observations.

%% Demo Parameters
N = 2^12;				% Signal dimensionality
subrate  = 0.60;		% Ratio of M/N (percent of dim. reduction)
sparsity = 0.50;		% Percent of signal which is non-zero
gb_mean = 0.2;			% Mean of GB signal prior
gb_var  = 0.75;			% Variance of GB signal prior
delta   = 1e-8;			% iid AWGN variance 
M = round(N*subrate);	% Number of measurements
K = round(sparsity*N);	% Number of non-zeros
xrange = 4*gb_var + gb_mean;

%% Generate Problem
A = randn(M,N) ./ sqrt(N);                    % A random iid projector		
x  = sqrt(gb_var).*randn(N,1) + gb_mean;	  % Generate gaussian part of signal...
rp = randperm(N);                             % Get random nonzero locations...
z = rp(K+1:end); 
nz = rp(1:K);
x(z) = 0;                                     % Set the zeros to make the signal sparse.
y = A*x + sqrt(delta)*randn(M,1);             % Calculate noisy measurements

%% Solve with ample-GB
fprintf('Running ample-GB...\n');
[a_gb,c_gb,history_gb] = ample(	A,y,@prior_gb,...
					     	   'prior_params', [0 1 0.5],...
                               'learn_prior_params',1,...
					  		   'true_solution',  x,...
                               'debug',0,...
                               'learn_delta',1, ...
                               'delta',1, ...
                               'convergence_tolerance',1e-10,...
                               'learning_mode','em',...
                               'report_history',1);

%% Reporting
fprintf('-----------------------\n');
mse_gb = norm(a_gb - x).^2./N;
fprintf('(GB) Final MSE : %0.2e\n',mse_gb);

%% Recovery Comparison
loc = 1:1:N;
figure(1); clf;
% Plot the ample-GB results
	hold on;
		stem(loc(nz),x(nz),'-bo','DisplayName','Original');
		stem(loc,a_gb,'-rx','MarkerSize',3,'DisplayName','Recovered Means');
	hold off;
	grid on; box on;
	xlabel(sprintf('MSE = %0.2e',mse_gb));
	title('Recovery with true GB prior');
	legend('Location','EastOutside');
	axis([1 N -xrange xrange]);

%% MSE Evolution Comparison
figure(2); clf;
        plot(history_gb.mse,'-b','LineWidth',2,'DisplayName','ample-GB');
    grid on; box on; 
    set(gca,'YScale','log');
    axis tight;
    title('MSE Evolution');
    legend('Location','NorthEast');

%% Delta Estimate Evolution Comparison
figure(3); clf;
    iters = length(history_gb.delta_estimate);
    dom = 1:1:iters;
    hold on;
        plot(history_gb.delta_estimate,'-b','LineWidth',2,'DisplayName','ample-GB');
        plot(dom,delta.*ones(iters,1),'-.k','LineWidth',1,'DisplayName','\Delta^*');
    hold off;
    grid on; box on; 
    set(gca,'YScale','log');
    axis([1 iters 1e-10 10]);
    title('Evolution of \Delta Estimation');
    legend('Location','NorthEast');

%% Prior Estimation Comparison
figure(4); clf;
    iters = size(history_gb.prior_params,1);
    dom = 1:1:iters;
    hold on;
        plot(history_gb.prior_params(:,1),'-r','LineWidth',2,'DisplayName','\mu');
        plot(history_gb.prior_params(:,2),'-b','LineWidth',2,'DisplayName','v');
        plot(history_gb.prior_params(:,3),'-g','LineWidth',2,'DisplayName','\rho');
        plot(dom,gb_mean.*ones(iters,1),'-k','LineWidth',1,'DisplayName','\mu^*');
        plot(dom,gb_var.*ones(iters,1),'-.k','LineWidth',1,'DisplayName','v^*');
        plot(dom,sparsity.*ones(iters,1),'--k','LineWidth',1,'DisplayName','\rho^*');
    hold off;
    grid on; box;
    axis tight;
    title('Evolution of Prior Parameters');
    legend('Location','EastOutside');