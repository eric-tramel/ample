% ample_qary_demo.m
%
% A demonstration test-case for ample.m. Will
% calculate the factorization of the posterior 
% distribution of an Q-valued 
% signal sampled with an iid Gaussian random matrix 
% given the set of observations. All learning will be
% done in EM mode, that is, update of learned parameters
% will be done after the AMP iteartion converges.

%% Demo Parameters
N = 2^12;				                % Signal dimensionality
subrate  = 0.30;		            % Ratio of M/N (percent of dim. reduction)
alphabet = 1:6;             % Alphabet of values
Na       = length(alphabet);
% % -- A well defined PMF for Na = 3.
% pmf       = [0.75,0.10,0.15];
% -- Control the "balance" of a more complicated PMF.
% -- One can observe how the more uniform, the harder it is for 
% -- ample to solve the problem. In this case, as maj_prob -> 1./Na;
maj_prob = 0.9;
pmf      = [maj_prob, (1-maj_prob)./(Na - 1) .* ones(1,Na-1) ];       % PMF of values
delta   = 1e-8;			            % iid AWGN variance 
em_iter = 40;                   % Max number of EM iterations
amp_iter = 1000;                % Max number of AMP iterations
M = round(N*subrate);	          % Number of measurements

%% Generate Problem
A = randn(M,N) ./ sqrt(N);                        % A random iid projector
x = alphabet(1).*ones(N,1);                       % Starting with a base signal
rp = randperm(N);                                 % Get random locations
start_idx = 1;
for i=1:length(alphabet)                          % Loop over the alpahbet...
  nAtValue = round(pmf(i)*N);                     % Determine the number of values according to pmf
  end_idx  = min(start_idx + nAtValue - 1,N);     % Make sure not to overrun the vector
  x(rp(start_idx:end_idx)) = alphabet(i);         % Set these values
  start_idx = end_idx + 1;
end
y = A*x + sqrt(delta)*randn(M,1);             % Calculate noisy measurements

%% Solve with ample-qary
fprintf('Running ample-qary...\n');
params{1} = alphabet;
params{2} = pmf;
[a_q,c_q,history_q] = ample(	A,y,@prior_qary,...
                               'prior_params', params,...
                               'learn_prior_params',0,...
                               'true_solution',  x,...
                               'debug',0,...
                               'learn_delta',1, ...
                               'delta',1, ...
                               'convergence_tolerance',1e-10,...
                               'learning_mode','track',...
                               'max_em_iterations',em_iter,...
                               'max_iterations',amp_iter, ...
                               'report_history',1,...
                               'convergence_type','iteration',...
                               'damp',0);

%% Reporting
fprintf('-----------------------\n');
mse_gb = norm(a_q - x).^2./N;
fprintf('  Final MSE : %0.2e\n',mse_gb);

%% Recovery Comparison
loc = 1:1:N;
figure(1); clf;
% Plot the ample-GB results
	hold on;
		stem(loc,x,'-bo','DisplayName','Original');
		stem(loc,a_q,'-rx','MarkerSize',3,'DisplayName','Recovered Means');
	hold off;
	grid on; box on;
	xlabel(sprintf('MSE = %0.2e',mse_gb));
	title('Recovery with Q-ary prior');
	legend('Location','EastOutside');
	axis([1 N 1-min(alphabet) 1+max(alphabet)]);

%% MSE Evolution Comparison
figure(2); clf;
        plot(history_q.mse,'-b','LineWidth',2,'DisplayName','ample-qary');
    grid on; box on; 
    set(gca,'YScale','log');
    axis tight;
    title('MSE Evolution');
    legend('Location','NorthEast');

%% Delta Estimate Evolution Comparison
figure(3); clf;
    iters = length(history_q.delta_estimate);
    dom = 1:1:iters;
    hold on;
        plot(history_q.delta_estimate,'-b','LineWidth',2,'DisplayName','ample-qary');
        plot(dom,delta.*ones(iters,1),'-.k','LineWidth',1,'DisplayName','\Delta^*');
    hold off;
    grid on; box on; 
    set(gca,'YScale','log');
    axis([1 iters 1e-10 10]);
    title('Evolution of \Delta Estimation');
    legend('Location','NorthEast');

% %% Prior Estimation Comparison
% figure(4); clf;
%     iters = size(history_q.prior_params,1);
%     dom = 1:1:iters;
%     hold on;
%         plot(history_q.prior_params,'-g','LineWidth',2,'DisplayName','\rho');
%         plot(dom,sparsity.*ones(iters,1),'--k','LineWidth',1,'DisplayName','\rho^*');
%     hold off;
%     grid on; box;
%     axis tight;
%     title('Evolution of Prior Parameters');
%     legend('Location','EastOutside');