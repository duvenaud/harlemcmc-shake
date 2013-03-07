function covs = rotate_cov( lengths, angle )
% Tamara Broderick
% David Duvenaud
%
% March 2013

unrotated = diag(lengths);
rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)];

covs = rotation' * unrotated * rotation;
