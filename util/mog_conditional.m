function cond_mix = mog_conditional( mix, x, fixed)
% Returns a mixture of Gaussians, conditioned on some of its
% dimensions being fixed to the value x.
%
% David Duvenaud
% March 2012

D = size(mix.means, 2);
notfixed = 1:D;
notfixed(fixed) = [];

cond_mix.weights = NaN(size(mix.weights, 1), 1);
cond_mix.means = NaN(size(mix.weights, 1), numel(notfixed));
cond_mix.covs = NaN(size(mix.weights, 1), numel(notfixed), numel(notfixed));

for k = 1:size(mix.means, 1);
    prec_all = inv(mix.covs(:,:,k));
    mu_a = mix.means(notfixed);
    mu_b = mix.means(fixed);
    prec_aa = prec_all(notfixed, notfixed, k);
    prec_ab = prec_all(notfixed, fixed, k);
    
    cond_mix.means(k) = mu_a - prec_aa \ prec_ab * ( x - mu_b);
    cond_mix.covs(k) = inv(prec_aa);
    evs = mvnpdf( x, mix.means(fixed), mix.covs(fixed, fixed, k));  % Evidences.
end

% Renormalize:
cond_mix.weights = mix.weights .* evs;
cond_mix.weights = cond_mix.weights ./ sum(cond_mix.weights);
