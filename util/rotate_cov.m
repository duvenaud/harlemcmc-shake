function covs = rotate_cov( lengths, angle )

unrotated = diag(lengths);
rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)];

covs = rotation' * unrotated * rotation;
