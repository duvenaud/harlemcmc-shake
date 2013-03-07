function [nll, dnll] = mix_gaussians_neg_log_pdf( x, mix )

[ll, dll] = mix_gaussians_log_pdf( x, mix );
nll = -ll;
dnll = -dll;



    

    

