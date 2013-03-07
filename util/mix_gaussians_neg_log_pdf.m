function [nll, dnll] = mix_gaussians_neg_log_pdf( x, mix )
% Tamara Broderick
% David Duvenaud
%
% March 2013

[ll, dll] = mix_gaussians_log_pdf( x, mix );
nll = -ll;
dnll = -dll;



    

    

