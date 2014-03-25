function [a,c,learned_params] = prior_l1sparse(r,s,params)
% PRIOR_L1SPARSE Generate the means and variances according to the
% 	prior parameters and the hidden variational variables {r,s} for
%	the assumption of the spare L1 prior (exponential distribution
% 	in the limit).
%
% [a,c] = prior_l1sparse(r,s,params) Calculate the means (a) and variances
%  (c) according to the given values of the prior parameters.
%  Params should consist of a 1x2 vector [min_x, max_x]. 
%     * If params is given as a Nx2 matrix, then it is assumed that the
%		signal is not iid and each value is calculated according to its 
%		own prior parameters.

%% Reassignments
x_min   = params(:,1);
x_max   = params(:,2);

%% Calculate means
a = max(x_min, (r + s) .* double((-r > s)) ) + ...
	min(x_max, (r - s) .* double(( r > s)) );

%% Calculate variances
c = s .* double(((r > s & r < (x_max + s)) | (r < -s & r > (x_min - s)) ));
c = max(c,1e-10);

%% Dummy Assignment
% There is no real change or "learning" of the parameters for this 
% prior, so we will just pass through the limits on x.
learned_params = params;