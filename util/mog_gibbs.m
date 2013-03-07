function sample = mog_gibbs( mix, x )

% Returns a sample from a mixture of Gaussians, conditioned on one of its
% dimensions being fixed.
%
% David Duvenaud
% Tamara Broderick
%
% March 2012

D = numel(x);

% Choose a dimension to sample.
dim = randi(D);
fixed = 1:D;
fixed(dim) = [];

% Find the conditional pdf.
cond_mix = mog_conditional( mix, x(fixed), fixed);

% Now sample from this distribution:
x(dim) = mix_gaussians_draw( cond_mix, 1 );

sample = x;
