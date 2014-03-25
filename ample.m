function [a,c,history] = ample(F,y,moment_func,varargin)
    [M,N] = size(F);
    F2 = F.^2;
    
    options = process_varargin(defaults(N,M),varargin);

    % Check to see if we need to learn parameters
    learn_prior = options.learn_prior_params;
    prior_params = options.prior_params;
    damp = options.damp;
    max_iter = options.max_iterations;
    conv_tol = options.convergence_tolerance;    
    delta = options.delta;    
    
    report_history = 0;
    if nargout > 2
        report_history = 1;
    end
    
    calc_mse = 0;
    if ~isempty(options.true_solution)
        calc_mse = 1;
        x0 = options.true_solution;
    end
    
    % Set initial loop parameters
    a = options.init_a;
    c = options.init_c;
    R = options.init_r;
    S = options.init_s;
    O = options.init_o;
    V = options.init_v;    

    %% Main Loop
    for i=1:max_iter
        % Update V
        oldDV = delta + V;
        V = F2*c;
        DV = delta + V;
        
        % Update \omega
        O = F*a - (y - O) .* (V./oldDV);

        % Update Sigma
        S_ = 1./(F2' * (1./DV));
        S = damp.*S + (1-damp).*S_;
        
        % Update R
        R_ = a + S.*(F'*((y - O)./DV));
        R = damp.*R + (1-damp).*R_;
        
        % Update moments
        last_a = a;
        if learn_prior
            [a,c,prior_params] = moment_func(R,S,prior_params);  
        else
            [a,c] = moment_func(R,S,prior_params);  
        end        
        
        % Update delta
        if options.learn_delta
            delta = sum((y - F*a).^2) ./ sum(1. ./ (1. + (F2*c) / delta));
        end
        convergence = norm(last_a - a).^2./N;
        
        % If we are in debug mode, show the current state
        if options.debug
            display_state(R,S,a,c);
        end

        % History Reporting
        if report_history
            history.convergence(i) = convergence;
            history.delta_estimate(i) = delta;
            history.prior_params(i,:) = prior_params;
            if calc_mse
                history.mse(i) = norm(a-x0).^2./N;
            end
        end

        % Check for convergence
        if convergence < conv_tol
            break;
        end
                
        % Output
        if calc_mse
            fprintf('\r [%d] | delta : %0.2e | convergence : %0.2e | mse : %0.2e      ',i,delta,convergence,norm(a-x0).^2./N);
        else
            fprintf('\r [%d] | delta : %0.2e | convergence : %0.2e |        ',i,delta,convergence);
        end
    end
    fprintf('\n');
    

function display_state(r,s,a,c)
    N = length(r);
    figure(2);
    subplot(2,2,1);
        stem(r,'-k.','MarkerSize',1);
        grid on; box on;
        axis tight;
        title('Varational Means');
    subplot(2,2,3);
        plot(s,'-k.','MarkerSize',1);
        grid on; box on;
        set(gca,'YScale','log');
        axis([1 N 1e-13 10]);
        title('Varational Variances');
    subplot(2,2,2);
        stem(a,'-k.','MarkerSize',1);
        grid on; box on;
        axis tight;
        title('Factorized Means');
    subplot(2,2,4);
        plot(c,'-k.','MarkerSize',1);
        grid on; box on;
        set(gca,'YScale','log');
        axis([1 N 1e-13 10]);
        title('Factorized Variances');
    drawnow;
        
    
function options = process_varargin(options,arguments)
    while ~isempty(arguments)
        field = lower(arguments{1});
        value = arguments{2};
        % Old Way
        %options = setfield(options,field,value);
        % New Way
        options.(field) = value;
        arguments(1:2) = [];
    end

function options = defaults(N,M)
    % AMPLE::DEFAULTS Set default options structure
    options.verbose_mode = 1;
    options.delta = 1;
    options.learn_delta = 1;
    options.prior_params = [];
    options.learn_prior_params = 0;
    options.max_iterations = 250;
    options.convergence_tolerance = 1e-10;
    options.damp = 0;
    options.true_solution = [];    
    options.init_v = ones(M,1);
    options.init_o = zeros(M,1);
    options.init_r = zeros(N,1);
    options.init_s = ones(N,1);
    options.init_a = zeros(N,1);
    options.init_c = ones(N,1);
    options.debug = 0;
    options.log_file = [];
