function [ll, dll] = mix_gaussians_log_pdf( x, mix )
% Evaluate a mixture of Guassians at specified locations.
%
% x is N x D
% mix.weights is k x 1
% mix.means is k x D
%
% David Duvenaud
% March 2012

% vals is going to be summed over the elements of mixture.
[K, D] = size(mix.means);
[N, D] = size(x);

log_pdfs = NaN(N, K);
log_mix_weights = log(mix.weights);
   
for k = 1:K
    log_pdfs(:, k) = log_mix_weights(k) + logmvnpdf( x, mix.means(k, :), mix.covs(:, :, k));
end
ll = logsumexp( log_pdfs );
    
if nargout > 1    
    d_log_pdfs = NaN(N, K, D);
    for k = 1:K
        [log_pdfs(:, k), d_log_pdfs(:, k, :)] = logmvnpdf( x, mix.means(k, :), mix.covs(:, :, k));
        d_log_pdfs(:, k, :) = d_log_pdfs(:, k, :) * mix.weights(k) * exp(log_pdfs(:, k));
    end
    
    dll = NaN( N, D);
    for d = 1:D
        dll(:, d) = sum( squeeze(d_log_pdfs( :, :, d)), 2) ./ exp(ll);
    end
end
