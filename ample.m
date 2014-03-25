function [a,c,history] = ample(F,y,moment_func,prior_params,x0)
    [M,N] = size(F);
    F2 = F.^2;

    damp = 0;
    max_iter = 300;
    conv_thresh = 1e-13;
    
    init_a = zeros(N,1);
    init_c = ones(N,1);  
    init_S = zeros(N,1);
    init_R = zeros(N,1);
    init_O = zeros(M,1);
    init_V = ones(M,1);
    delta = 1;    
    
    report_history = 0;
    if nargout > 2
        report_history = 1;
    end
    
    debug = 0;
    if nargin > 4
        debug = 1;
    end
    
    a = init_a;
    c = init_c;
    O = init_O;
    V = init_V;
    R = init_R;
    S = init_S;

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
        [a,c] = moment_func(R,S,prior_params);  
        
        % Update delta
        delta = sum((y - F*a).^2) ./ sum(1. ./ (1. + (F2*c) / delta));
        convergence = norm(last_a - a).^2./N;

        % History Reporting
        if report_history
            history.convergence(i) = convergence;
            history.delta_estimate(i) = delta;
            state.a = a;
            state.c = c;
            state.r = R;
            state.b = 1./S;
            state.D = delta.*ones(M,1);            
            if debug
                history.mse(i) = norm(a-x0).^2./N;
            end
        end

        % Check for convergence
        if convergence < conv_thresh
            break;
        end
                
        % Output
        if debug
            fprintf('\r [%d] | delta : %0.2e | convergence : %0.2e | mse : %0.2e      ',i,delta,convergence,norm(a-x0).^2./N);
        else
            fprintf('\r [%d] | delta : %0.2e | convergence : %0.2e |        ',i,delta,convergence);
        end
    end
    fprintf('\n');
    