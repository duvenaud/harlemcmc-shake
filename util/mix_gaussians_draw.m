function vals = mix_gaussians_draw( mix, N )
% Generate random draws from a mixture of Gaussians.
%
% Tamara Broderick
% David Duvenaud
%
% March 2013

vals = NaN(N, size(mix.means,2));

for i = 1:N;
    cur_component = find(mnrnd( 1, mix.weights ));
    vals(i, :) = mvnrnd( mix.means(cur_component, :), ...
                         mix.covs(:, :, cur_component));
end
