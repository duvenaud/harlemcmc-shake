function vals = mix_mix_gaussians_pdf( x, mixes )
% Evaluate a mixutre of mixtures of Guassians at specified locations.
%
% Tamara Broderick
% David Duvenaud
%
% March 2013

vals = zeros(size(x, 1), 1);
for n = 1:numel(mixes)
    mix = mixes{n}
    for k = 1:size(mix.means, 1);
        vals = vals + mix.weights(k) .* mvnpdf( x, mix.means(k, :), mix.covs(:, :, k));
    end
end
