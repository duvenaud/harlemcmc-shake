function vals = mix_gaussians_pdf( x, mix )
% Evaluate a mixture of Guassians at specified locations.
%
% Tamara Broderick
% David Duvenaud
%
% March 2013

vals = zeros(size(x, 1), 1);
for k = 1:size(mix.means, 1);
    vals = vals + mix.weights(k) .* mvnpdf( x, mix.means(k, :), mix.covs(:, :, k));
end
