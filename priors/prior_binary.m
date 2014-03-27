function [a,c,learned_params] = prior_binary(r,s,params)
% PRIOR_BINARY Generate the means and variances according to the
% 	prior parameters and the hidden variational variables {r,s}.
%
% [a,c,learned_params] = prior_binary(r,s,params) Calculate the means (a) and variances
%  (c) according to the given values of r, s, and the prior parameters.
%  Params should consist of a 1-value vector, [rho]. 
%     * If params is given as a Nx1 vector, then it is assumed that the
%		signal is not iid and each value is calculated according to its 
%		own prior parameters.
%     * The value s should be given as the square.

%% I/O
learn_prior = 0;
if nargout > 2
    learn_prior = 1;
end

%% Reassignments
n   = length(r);
rho = params;

%% Partition calculation
z = rho + (1-rho) .* exp(0.5 .* (1 - 2.*r)./s);
a = rho ./ z;
c = a .* (1 - a);
                                    
%% Learn New Prior 
if learn_prior
	if size(params,1) ~= 1
		error('Prior learning not implemented for non-iid signals.\n');
	end
	% Update Rho
	rho = mean(a);

	learned_params = [rho];
end
