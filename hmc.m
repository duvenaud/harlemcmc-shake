function [params, nll, arate, tail] = hmc(likefunc, x, options, varargin)     
% Hamiltonian Monte Carlo - copied from David Mackay's book.
%
% David Duvenaud
% Tomoharu Iwata
%
% April 2012
%
% likefunc returns nll, dnll
%
% options.Tau is the number of leapfrog steps. 
% options.epsilon is step length


silent = true;

arate = 0; %acceptance rate
L = options.num_iters;

[E, g] = likefunc( x, varargin{:});



for l = 1:L
    p = randn( size( x ) );
    H = p' * p / 2 + E;
    
    xnew = x; gnew = g;
    
    % Randomize step length and number of steps.
    %cur_tau = randi(options.Tau);
    %cur_eps = rand * options.epsilon;
    
    cur_tau = options.Tau;
    cur_eps = options.epsilon;
    
    tail = NaN(cur_tau+1, size(x,2));
    tail(1,:) = x;
    
    for tau = 1:cur_tau
        p = p - cur_eps * gnew / 2;
        xnew = xnew + cur_eps * p;
        tail(tau+1, :) = xnew;
        [ignore, gnew] = likefunc( xnew, varargin{:}); 
        
        p = p - cur_eps * gnew / 2;
    end
    
    [Enew, ignore] = likefunc( xnew, varargin{:});    
    Hnew = p' * p / 2 + Enew;
    dh = Hnew - H;
    
    if dh < 0
        accept = 1;
        if ~silent; fprintf('a'); end
    else
        if rand() < exp(-dh)
            accept = 1;
            if ~silent; fprintf('A'); end
        else
            accept = 0;
            if ~silent; fprintf('r'); end
        end
    end
    
    if accept
        g = gnew;
        x = xnew;
        E = Enew;
        arate = arate+1;
    else
        tail = [];
    end
end
 
arate = arate/L;
params = x;
nll = E;
