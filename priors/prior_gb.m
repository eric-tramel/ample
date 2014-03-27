function [a,c,learned_params] = prior_gb(r,s,params)
% PRIOR_GB Generate the means and variances according to the
% 	prior parameters and the hidden variational variables {r,s}.
%
% [a,c] = prior_gb(r,s,params) Calculate the means (a) and variances
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
n   = length(r);
m   = params(:,1);
v   = params(:,2);
rho = params(:,3);

%% Partition calculation
r2 = r.*r;
vps = v + s;
mmr = m - r;
mmr2 = mmr.*mmr;
fac = rho .* sqrt(s ./ vps);
rsc = -0.5 .* mmr2 ./ vps;    
zIn = ((1-rho)./fac) .* exp(-0.5.*r2./s - rsc);
z = 1 + zIn;

%% Moment Calculation  
crs = m.*s + r.*v;
crs2 = crs.*crs;
vpr = v.*s.*vps;
ivps = 1./ vps;
ivps2 = ivps.*ivps;

a = ivps .* crs ./ z;
a2 = a.*a;

m2 = ivps2 .* (crs2 + vpr)./z;
c = max(1e-18,m2 - a2);				% Ensure that we don't get a numerically bad
									% value for c.
                                    
%% Learn New Prior 
if learn_prior
	if size(params,1) ~= 1
		error('Prior learning not implemented for non-iid signals.\n');
	end
	% Update Rho
	cs = r ./ s + m ./ v;
	sv = 1 ./ s + 1 ./ v;
	ex = exp(0.5 .* cs .* cs ./ sv  -  0.5 .* m .* m ./ v);
	rho_n = sum(sv.*a./cs);
	rho_d = sum( 1 ./ ((1-rho) + rho .* ex ./ sqrt(v .* sv)) );
	rho = rho_n ./ rho_d;

	% Update Mean 
	m = sum(a) ./ (rho.*n);

	% Update Variance
	v = sum(c + a.*a) ./ (rho * n) - m.*m; 

	learned_params = [m, v, rho];
end
