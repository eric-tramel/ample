function [a,c,learned_params] = prior_complexgb(r,s,params)
% PRIOR_COMPLEXGB Generate the means and variances according to the
% 	prior parameters and the hidden variational variables {r,s}.
%
% [a,c] = prior_complexgb(r,s,params) Calculate the means (a) and variances
%  (c) according to the given values of r, s, and the prior parameters.
%  Params should consist of a 1x3 vector [mean,variance,rho]. 
%     * If params is given as a Nx3 matrix, then it is assumed that the
%		signal is not iid and each value is calculated according to its 
%		own prior parameters.
%     * The value s should be given as the square.

%% I/O
learn_prior = 0;
if nargout > 2
    learn_prior = 1;
end

%% Reassignments
m   = params(:,1);
v   = params(:,2);
rho = params(:,3);

%% Moment Calculation
% From Jean's Code
%             var_gauss = prior.param_2; mean_gauss = prior.param_1; R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho;
%             prior.av_mess_old = prior.av_mess;
%             chi2 = var_gauss .* S2_./ (var_gauss + S2_);
%             M = (var_gauss .* R_ + S2_ .* mean_gauss) ./ (S2_ + var_gauss);
%             alpha_ = abs(mean_gauss).^2 ./ var_gauss + abs(R_).^2 ./ S2_ - abs(M).^2 ./ chi2;
%             Z = (1 - rho_) .* exp(-0.5 .* abs(R_).^2 ./ S2_) + rho_ .* S2_ ./ (S2_ + var_gauss) .* exp(-alpha_ ./ 2);
%             prior.av_mess = (rho_ .* S2_ ./ (S2_ + var_gauss) .* M .* exp(-alpha_ ./ 2) ) ./ Z;
%             prior.var_mess = max(1e-18, (rho_ ./ Z .* S2_ ./ (S2_ + var_gauss) .* (2 .* chi2 + abs(M).^2) .* exp(-alpha_ ./ 2) - abs(prior.av_mess).^2) ./ 2);
%             prior.var_mess(~isfinite(prior.var_mess) ) = 0;
chi2   = v .* s ./ (m + s);
M      = (v.*r + s.*m) ./ (s + m);
alpha_ = abs(m).^2 ./ v + abs(r).^2./s - abs(M).^2 ./ chi2;
Z      = (1-rho) .* exp(-0.5 .* abs(r).^2./s) + rho.* s ./ (s+v) .*exp(-alpha_ ./ 2);
a      = (rho .* s ./ (s+v) .* M .* exp(-alpha_./2))./Z;
c      = (rho ./Z .* s./(s+v) .* (2 .* chi2 + abs(M).^2) .* exp(-alpha_ ./ 2) - abs(a).^2)./2;
c(~isfinite(c)) = 0;
c      = max(c,1e-18);

%% Learn New Prior 
if learn_prior
	if size(params,1) ~= 1
		error('Prior learning not implemented for non-iid signals.\n');
    end
    % The learning of the prior in this complex case has not yet been
    % calculated.
    learned_params = params;
end
