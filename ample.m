function [a,c,history,R,S] = ample(F_,y,moment_func,varargin)
    mean_approximation = 0;
    if ~isa(F_,'struct')
    % In this case, the user has supplied us with an explicit matrix for
    % the system we are trying to solve.
        [M,N] = size(F_);
        % Storing these in memory to save on computation.        
        FT_  = F_';
        F2_  = abs(F_).^2;
        F2T_ = F2_';
        % Assignments to the anonymous functions.
        F   = @(x_) F_*x_;
        FT  = @(x_) FT_*x_;
        F2  = @(x_) F2_*x_;
        F2T = @(x_) F2T_*x_;        
    else
    % In this case, the user has supplied us with a structure continaing
    % the forward, squared-forward, adjoint, and squared-adjoint operators.
    % This structure must additionally contain the signal dimensionality
    % and number of measurements.
        N = F_.N;
        M = F_.M;        
        
        % Assignments
        F   = F_.forward;
        FT  = F_.adjoint;
        % Sometimes we may not actually have these squared forms, so we
        % will only set these if they exist.
        if ~isfield(F_,'squared_forward') && ~isfield(F_,'squared_adjoint')
            % If the user hasn't specified a squared-forward and -adjoint
            % operator, we can only run in the unity-approx mode...            
            mean_approximation = 1;
            % ... but we should probably warn them about it.
            fprintf('[WARNING] F.squared_forward or F.squared_adjoin unspecified. Running in unity-approximation mode.\n');
        end        
        if isfield(F_,'squared_forward')
            F2  = F_.squared_forward;
        end
        if isfield(F_,'squared_adjoint')            
            F2T = F_.squared_adjoint;
        end
    end
    
    % Set the default options and update them according to the
    % user specified options.
    options = process_varargin(defaults(N,M),varargin);

    % Re-assignments to put off having to refactor the code.
    mean_approximation = mean_approximation | options.mean_approximation;    
    learn_prior = options.learn_prior_params;
    prior_params = options.prior_params;
    damp = options.damp;
    max_iter = options.max_iterations;
    conv_tol = options.convergence_tolerance;    
    delta = options.delta;    

    em_mode = 1;
    switch(options.learning_mode)
        case 'em'
            em_mode = 1;
        case 'track'
            em_mode = 0;
    end
    
    report_history = 0;
    history = [];
    if nargout > 2 && options.report_history
        report_history = 1;
    end
    
    calc_mse = 0;
    if ~isempty(options.true_solution)
        calc_mse = 1;
        x0 = options.true_solution;
    end
    
    if mean_approximation
        F2  = @(x_) mean(x_).*ones(M,1);
        F2T = @(x_) mean(x_).*ones(N,1);
    end
    
    % Set initial loop parameters
    a = options.init_a;
    c = options.init_c;
    R = options.init_r;
    S = options.init_s;
    O = options.init_o;
    V = options.init_v;    

    if ~em_mode
    % In this case, we solve the problem using parallel learning updates, or
    % no learning updates at all. This is the default operation of ample.
        for i=1:max_iter
            % Update V
            oldDV = delta + V;
            V = F2(c);
            DV = delta + V;
            
            % Update \omega
            O = F(a) - (y - O) .* (V./oldDV);

            % Update Sigma
            S_ = 1./F2T(1./DV);
            S = damp.*S + (1-damp).*S_;
            
            % Update R
            R_ = a + S.*(FT((y - O)./DV));
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
                delta = sum( (abs(y) - abs(O)).^2  ./ (1 + V./delta).^2 ) ./ sum(1 ./ (1 + V./delta));
            end
            
            convergence = norm(last_a - a).^2./N;
                    
            % History Reporting
            if report_history
                history.convergence(i) = convergence;
                history.delta_estimate(i) = delta;
                history.prior_params(i,:) = prior_params;
                if calc_mse
                    history.mse(i) = norm(a-x0).^2./N;
                end
            end

            % If we are in debug mode, show the current state
            if options.debug
                display_state(R,S,a,c);
            end

            % Do we have nan issues?
            nancheck(V,'V','error');
            nancheck(O,'O','error');
            nancheck(R,'R','error');
            nancheck(S,'S','error');
            nancheck(a,'a','error');
            nancheck(c,'c','error');


            % Output
            if calc_mse
                fprintf('\r [%d] | delta : %0.2e | convergence : %0.2e | mse : %0.2e      ',i,delta,convergence,norm(a-x0).^2./N);
            else
                fprintf('\r [%d] | delta : %0.2e | convergence : %0.2e |        ',i,delta,convergence);
            end
            
            % Check for convergence
            if convergence < conv_tol
                break;
            end                
        end
    else
    % In this case, we want to use EM-learning, which means that we update the learned parameters
    % only after convegence of the AMP iteartion. We can accomplish this by making a recursive
    % call to ample.
        emoptions = options;
        emoptions.learn_delta = 0;                  % Turn off learning mode
        emoptions.learn_prior_params = 0;           % Turn off learning mode
        emoptions.learning_mode = 'track';          % Make sure we don't EM inside the EM
        % Need to initialize the history variables
        if report_history
            history.convergence = [];            
            history.delta_estimate = [];
            history.prior_params = [];
            if calc_mse
                history.mse = [];
            end
        end
        for emiter = 1:options.max_em_iterations
            fprintf('**EM ITERATION #%d**\n',emiter);
            % Call ample for this local solution
            last_a = a;
            input_args = struct2varargin(emoptions);
            [a,c,history_,r,s] = ample(F_,y,moment_func,input_args{:});

            % Update the prior parameters
            if options.learn_prior_params
                % Need to have the final R and S values...
                [~,~,emoptions.prior_params] = moment_func(r,s,emoptions.prior_params);  
            end

            % Update Delta -- Using the fixed-point omega style update
            if options.learn_delta
                emoptions.delta = sum(abs(y - F(a)).^2) ./ sum(1 ./ (1. + F2(c)./emoptions.delta));
            end

            % TODO: Merging of the per-EM-iteration history terms.
            if report_history
                history.convergence = [history.convergence, history_.convergence];
                history.delta_estimate = [history.delta_estimate, history_.delta_estimate];
                history.prior_params = [history.prior_params; history_.prior_params];
                if calc_mse
                    history.mse = [history.mse, history_.mse];
                end
            end
            
            % Update the initial states
            emoptions.init_a = a;
            emoptions.init_c = c;
            emoptions.init_r = r;
            emoptions.init_s = s;

            % Check for convergence
            convergence = norm(a - last_a).^2./N;
            if convergence < conv_tol
                break;
            end            
        end
    end
    fprintf('\n');
    
%% Helper Functions
function display_state(r,s,a,c)
    % AMPLE::DISPLAY_STATE Visualization function when running in debug
    % mode. Displays the current values of {r,s,a,c}. Plots the
    % variance-type variables, {s,c}, on a log-scale.
    N = length(r);
    figure(42);
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
    % AMPLE::PROCESS_VARARGIN Scan through all the user arguments and set
    % the fields of the options structure accordingly. Assumes that the
    % names the user sets and the target fields in the structure are the
    % same.
    while ~isempty(arguments)
        field = lower(arguments{1});
        value = arguments{2};
        options.(field) = value;
        % If having compatibility issues, comment out the above line and
        % use the subsequent one, instead.
        %options = setfield(options,field,value);               
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
    options.mean_approximation = 0;
    options.learning_mode = 'track';
    options.max_em_iterations = 20;
    options.report_history = 1;
    
function vars = struct2varargin(structure)
    % AMPLE::STRUCT2VARARGIN Convert the given structure to varargin
    % format, that is, all odd-indexed cells containing the field names and
    % all even-indexed cells containing the values.
    var_names = fieldnames(structure);
    var_vals  = struct2cell(structure);
    vars = cell(2*length(var_names),1);
    vars(1:2:end) = var_names;
    vars(2:2:end) = var_vals;

function nancheck(x,name,report_level)
    % AMPLE::NANCHECK Checks the value of x. If it is a NaN, an report
    % will be made.
    nanid = isnan(x);
    snid = sum(nanid);
    if snid > 0
        n = length(x);
        rep_str = sprintf('NaN dectected on %d entries of (%s), [size=%d].',snid,name,n);
        switch report_level
            case 'warning'
                fprintf('[AMPLE::WARNING] %s \n',rep_str);
            case 'error'
                error('[AMPLE::ERROR] %s \n',rep_str);
        end
    end
    
    
    