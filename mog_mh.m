function sample = mog_mh( mix, cur_pt, proposal_cov )
%
% Metropolis-Hastings sampling in a mixture of Gaussians.
%
% Tamara Broderick
% David Duvenaud
% March 2013

% Compute MH proposal.
proposal = mvnrnd( cur_pt, proposal_cov );
proposal_ll = mix_gaussians_log_pdf(proposal, mix);
cur_ll = mix_gaussians_log_pdf(cur_pt, mix);

% Possibly take a MH step.
ratio = exp(proposal_ll - cur_ll);
if ratio > rand
    sample = proposal;   % Accept. :)
else
    sample = cur_pt;     % Reject. :(
end
