function [a,c,learned_params] = prior_qary(r,s,params)
% PRIOR_QARY Generate the means and variances according to the
% 	prior parameters and the hidden variational variables {r,s}.
%
% [a,c,learned_params] = prior_qary(r,s,params) Calculate the means (a) and variances
%  (c) according to the given values of r, s, and the prior parameters.
%  Params should consist of a length-q vector of the PMF of each value.
%     * The value s should be given as the square.
%	  * The prior has the following format...
%	    params{1}   :   Alphabet of possible values (Qx1)
%		params{2}	: 	PMF of value abundances, in same order as params{1} (Qx1 vector.)
%						Must sum to 1.

%% I/O
learn_prior = 0;
if nargout > 2
    learn_prior = 1;
end

fprintf('max R = %0.2e',max(r));
%% Reassignments
n   = length(r);
alphabet = params{1};
pmf      = params{2};
nv   = length(pmf);
exp_part = zeros(n,nv);
z  = zeros(n,1);
a  = zeros(n,1);
a2 = zeros(n,1);

%% Ratio Approach
% 1. Find out the maximum un-normalized "probability"
p_tilde_score = zeros(nv,n);
for i=1:nv
	p   = pmf(i);
	tau = alphabet(i);
	p_tilde_score(i,:) = (2.*(tau.*r) - tau*tau)./(2.*s) + log(p);
end
[~, jmax] = max(p_tilde_score);
taumax  = alphabet(jmax);
taumax2 = taumax.*taumax;
taumax = taumax(:);
taumax2 = taumax2(:);
p0max   = pmf(jmax);

% 2. Calculate the "partition" according to the ratio against pmax...
for i=1:nv
	p    = pmf(i);
	tau  = alphabet(i);
	tau2 = tau.*tau;
	logpr = log(p./p0max);
    logpr = logpr(:);
	pratio = exp((2.*tau.*r - 2.*taumax.*r + taumax2 - tau2)./(2.*s) + logpr); 
	z = z + pratio;	
end

% 3. Calculate "a" and "a2"
for i=1:nv
	p    = pmf(i);
	tau  = alphabet(i);
	tau2 = tau.*tau;
	logpr = log(p./p0max);
    logpr = logpr(:);
	pratio = exp((2.*tau.*r - 2.*taumax.*r + taumax2 - tau2)./(2.*s) + logpr); 
	a  = a  + (tau .*pratio)./z;
	a2 = a2 + (tau2.*pratio)./z;	
end

% 4. Calculate the variances, c
c = max(a2 - a.*a,1e-18);

                                    
%% Learn New Prior 
if learn_prior
	if size(params,1) ~= 1
		error('Prior learning not implemented for non-iid signals.\n');
	end

	fprintf('\n[WARNING] Prior learning is not completely verified for Q-ary.\n');

	% TODO: Learn the pmf...
	% 1. Calculate the MAP estimate according to the {R,S} values.
	q = zeros(nv,n);
	for i=1:nv
		p    = pmf(i);
		tau  = alphabet(i);
		tau2 = tau.*tau;
		logpr = log(p./p0max);
		logpr = logpr(:);
		pratio = exp((2.*tau.*r - 2.*taumax.*r + taumax2 - tau2)./(2.*s) + logpr); 
		q(i,:) = pratio ./ z;
	end
	[~,max_states] = max(q);
	x = alphabet(max_states);

	% 2. Count the ratio of each dictionary value in the MAP estimate
	%    to determine the empirical PMF.
	x_alphabet = unique(x);
	for i=1:length(x_alphabet)
		new_pmf(i) = sum(x == x_alphabet(i)) ./ n;
	end
	learned_params     = params;
	learned_params{2}  = new_pmf;
end