function harlem_shake( save_pngs_flag )
% Harlem shake - MCMC version.
%
% Tamara Broderick
% David Duvenaud
%
% March 2013


addpath('exportfig');
addpath('util');


mh_time = 16;   % Seconds of buildup part.
hmc_time = 23;  % Seconds of crazy part.
framerate = 14; % Hz.
n_mh_frames = mh_time * framerate;
n_hmc_frames = hmc_time * framerate;

% Check for previously cached samples.
cache_filename = 'sample_cache.mat';
if exist(cache_filename, 'file')
    load(cache_filename);
else
    mh_proposal_cov = [ 0.005 0; 0 0.005 ];
    mh_proposal_middle_cov = [ 0.5 0; 0 0.5 ];

    % Set up some mixtures of Gaussians.
    mixes = define_mixes();
    num_mixes = numel(mixes);
    
    % Fix the seed of the random generators.
    seed=0;
    randn('state',seed);
    rand('state',seed);
    
    % Start with samples from these dists.
    x = cell(num_mixes, 1);
    for n = 1:num_mixes
        x{n} = mix_gaussians_draw( mixes{n}, 1 );
    end    
    
    % Run M-H.
    fprintf('\nComputing MH samples');
    for m = 1:num_mixes
        samples{m} = NaN(n_mh_frames, 2);
        samples{m}(1,:) = x{m};
        
        if m == 5   % The central Gaussian gets its own proposal distribution.
            cur_cov = mh_proposal_middle_cov;
        else
            cur_cov = mh_proposal_cov;
        end
        for n = 2:n_mh_frames
            samples{m}(n,:) = mog_mh( mixes{m}, samples{m}(n - 1,:), cur_cov );
        end
        fprintf('.');
    end

    % Fix the seed of the random generators.
    randn('state',seed);
    rand('state',seed);

    fprintf('\nComputing HMC samples');
    for m = 1:num_mixes

        if 0
            % NUTS.
            Madapt = 250;
            loglikefunc = @(t) mix_gaussians_log_pdf( t, mixes{m} );
            lambda = 0.1;
            [samples{m}, epsilon] = hmc_da(loglikefunc, n_frames, Madapt, x{m}, lambda);
        else
            % Mackay brand HMC
            hmc_options.num_iters = 1;
            hmc_options.Tau = 20;  % Number of steps.
            hmc_options.epsilon = 0.05;
            hmc_samples{m} = NaN(n_hmc_frames, 2);
            hmc_samples{m}(1,:) = x{m};
            for n = 2:n_hmc_frames
                loglikefunc = @(t) mix_gaussians_neg_log_pdf( t, mixes{m} );
                [hmc_samples{m}(n,:), nll, arate, tail{m,n}] = hmc( loglikefunc, hmc_samples{m}(n - 1,:), hmc_options );
            end
        end

        fprintf('.');
    end
    save(cache_filename);
end

if nargin < 1
    save_pngs = false; 
else
    save_pngs = save_pngs_flag;
end

visual_framerate = framerate;
plot_hmc_tails = true;


figure(1); clf;
set(gcf, 'Position',[1 1 1400 1050]);

% Plot MH part.
frame_number = 1;
frame_number = plot_samples( samples, n_mh_frames, num_mixes, save_pngs, visual_framerate, frame_number, mixes, [], plot_hmc_tails );

fprintf('\n\nDO THE HARLEM SHAKE\n\n')

% Plot HMC part.
plot_samples( hmc_samples, n_hmc_frames, num_mixes, save_pngs, visual_framerate, frame_number, mixes, tail, plot_hmc_tails );


% Compile the video.
if save_pngs
    system('ffmpeg -r 14 -i frames/hs_%04d.png -vcodec huffyuv hs_movie_v7.avi');
end

end


function frame_number = plot_samples( samples, n_frames, num_mixes, ...
                             save_pngs, framerate, frame_number, mixes, tail, plot_hmc_tails )

history = 10;
num_hmc_tails = 9;
col_change_frames = 6;


in_hmc = numel(tail) > 0;
if in_hmc
    ls = 'none';
else
    ls = '-';
end

% Cache contours.
for m = 1:num_mixes
    margin = 0.01;
    h_axes(m) = subaxis(3,3,m,'Spacing',0.01, 'MR',0.01, 'Holdaxis', true, ...
        'MarginLeft',margin,'MarginRight',margin, ...
        'MarginTop',margin,'MarginBottom',margin);
    
    plot_one_contour(mixes{m}); hold on;
    set(gca, 'LooseInset', [0,0,0,0]);
end

col_ix = 2;
c_array(1, :) = [ 55, 126, 184 ];  % blue
c_array(2, :) = [ 255, 127, 1 ];   % orange
c_array(3, :) = [ 77, 175, 74 ];   % green
c_array(4, :) = [ 250, 60, 80 ];   % red
c_array(5, :) = [ 152, 78, 163 ];  % purple
c_array(6, :) = [ 200, 255, 51 ];  % yellow


for n = 1: n_frames
    
    % Change the background color if we're in the HMC phase.
    if in_hmc && (mod(n, col_change_frames) == 1)
        set(gcf, 'color', c_array(mod(col_ix, 6) + 1, :) ./255);
        col_ix = col_ix + 1;
    end

    cur_range = max(1, n - history):n;
    cur_tail_range = max(1, n - num_hmc_tails):n;
    
    % Plot Gibbs samplers running.
    for m = 1:num_mixes
        % Plot the sample.
        if plot_hmc_tails && in_hmc
            % Show HMC tail.
            for t_ix = cur_tail_range
                if numel(tail{m,t_ix}) > 0
                    h_tail{m, t_ix} =  plot(h_axes(m), tail{m,t_ix}(:,1), tail{m,t_ix}(:,2), ...
                        'c-', 'LineWidth', 5, 'Color', colorbrew(6));
                end
            end
        end
        
        h{m} = plot(h_axes(m), samples{m}(cur_range,1), samples{m}(cur_range,2), ...
            '-', 'LineWidth', 7, 'Marker', 'o', 'MarkerSize', 15, ...
            'LineStyle', ls, 'MarkerFaceColor', 'r', 'Color', colorbrew(6), ...
            'MarkerEdgeColor', 'r');
    end
    
    pause(1/framerate);

    if save_pngs
        set(gcf, 'Position',[1 1 1024 768]);
        export_fig('-nocrop', sprintf('frames/hs_%04d.png', frame_number));
    end
    frame_number = frame_number + 1;
    
    % Erase old dots.
    for m = 1:num_mixes
        delete(h{m});
        if plot_hmc_tails && in_hmc
            for t_ix = cur_tail_range
                if numel(tail{m,t_ix}) > 0
                    delete(h_tail{m, t_ix});
                end
            end
        end        
    end
end

end



function dist_vals = plot_one_contour(mix, dist_vals)
    % Plot the contours.
    length = 2;
    range = [ -length, length; -length length];
    N_1d = 100;
    ncontours = 4;

    xrange = linspace( range(1,1), range(1,2), N_1d);   % Choose a set of x locations.
    yrange = linspace( range(2,1), range(2,2), N_1d);   % Choose a set of x locations.
    [xvals, yvals] = meshgrid( xrange, yrange);
    gridvals = [xvals(:) yvals(:)];
       
    if nargin < 2
        dist_vals = mix_gaussians_pdf(gridvals, mix );
    end
    
    %colormap('Gray');
    %map = [0 0 0; 1 1 1];
    %sc = 0.5; ec = 1;
    %gradient = linspace( sc, ec, 100)'.^1;
    %map = [gradient, gradient, gradient];
    map = [1 1 1];
    colormap(map);
    
    dh = contour( xvals, yvals, reshape(dist_vals .^ 0.6, N_1d, N_1d ), ncontours, ...
        'LineWidth', 1); hold on;

    % Make plot prettier.
    set(gcf, 'color', 'white');
    set(gca, 'color', 'black');
    set(gca, 'YGrid', 'off');
    set(gca, 'Xtick', []);
    set(gca, 'Ytick', []);
    %axis off
end



function mixes = define_mixes()

mixes = cell(0);

% H.
if 1
skinny = 0.01;
really_skinny = 0.005;
fat = 0.6;
vert_cov = [really_skinny 0; 0 fat];
horz_cov = [fat 0; 0 skinny];
mix.means = [ -1.5 0; 0 0; 1.5 0];
mix.covs(:,:,1) = vert_cov;
mix.covs(:,:,2) = horz_cov;
mix.covs(:,:,3) = vert_cov;
mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mixes{end + 1} = mix;
end

% A.
if 1
skinny = 0.01;
fat = 0.45;
horz_cov = [fat 0; 0 skinny];
mix.means = [ 0 -1; -0.9 0; 0.9 0];
mix.covs(:,:,1) = horz_cov;
lengths = [0.01 1.2];
mix.covs(:,:,2) = rotate_cov( lengths, pi/8 );
mix.covs(:,:,3) = rotate_cov( lengths, -pi/8 );
mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mixes{end + 1} = mix;
end


% R.
if 1
skinny = 0.01;
fat = 0.75;
vert_cov = [skinny 0; 0 fat];
mix.means = [ -1.5 0; 0.25 1.25; 0.25 0.5; 0.25 -0.9 ];
mix.covs(:,:,1) = vert_cov;
lengths = [1.2 0.01];
mix.covs(:,:,2) = rotate_cov( lengths, pi/12 );
mix.covs(:,:,3) = rotate_cov( lengths, -pi/12 );
mix.covs(:,:,4) = rotate_cov( lengths, pi/8 );
mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mixes{end + 1} = mix;
end

% L.
if 1
skinny = 0.01;
fat = 0.6;
vert_cov = [skinny 0; 0 fat];
horz_cov = [fat 0; 0 skinny];
mix.means = [ -1.5 0; 0 -1.5];
mix.covs(:,:,1) = vert_cov;
mix.covs(:,:,2) = horz_cov;
mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mixes{end + 1} = mix;
end

% Spherical Gaussian.
mix.weights = 1;
mix.means = [ 0 0];
mix.covs = [ 1 0; 0 1];
mixes{end + 1} = mix;


% M.
if 1
skinny = 0.01;
really_skinny = 0.005;
fat = 0.6;
vert_cov = [really_skinny 0; 0 fat];
mix.means = [ -1.5 0; 1.5 0; -0.7 0.2; 0.7 0.2 ];
mix.covs(:,:,1) = vert_cov;
mix.covs(:,:,2) = vert_cov;
lengths = [0.01 1];
mix.covs(:,:,3) = rotate_cov( lengths, -pi/6 );
mix.covs(:,:,4) = rotate_cov( lengths, pi/6 );
mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mixes{end + 1} = mix;
end

%S
if 1
num_circle = 50;
angles = linspace( -pi/2, pi, num_circle) + pi/2;
mix.means = [cos( angles').*1.75 - 1, sin(angles').*0.9] + 0.9;
mix.means = [mix.means; -mix.means];
mix.covs = repmat( [ 1 0; 0 1], [1, 1, num_circle*2]) .* 0.01;
mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mixes{end + 1} = mix;
end

% H.
if 1
skinny = 0.1;
really_skinny = 0.005;
fat1 = 1.2;
fat2 = 0.6;
vert_cov = [really_skinny 0; 0 fat1];
horz_cov = [fat2 0; 0 skinny];
mix.means = [ -1.5 0; 0 0; 1.5 0];
mix.covs(:,:,1) = vert_cov;
mix.covs(:,:,2) = horz_cov;
mix.covs(:,:,3) = vert_cov;
%mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mix.weights = [0.1 0.8 0.1];
mixes{end + 1} = mix;
end


% K.
if 1
skinny = 0.01;
fat = 0.75;
vert_cov = [skinny 0; 0 fat];
mix.means = [ -1.5 0; -0.1 -0.8; -0.1 0.8 ];
mix.covs(:,:,1) = vert_cov;
lengths = [1.2 0.01];
mix.covs(:,:,2) = rotate_cov( lengths, pi/6 );
mix.covs(:,:,3) = rotate_cov( lengths, -pi/6 );
mix.weights = ones(size(mix.means,1),1) ./ size(mix.means,1);
mixes{end + 1} = mix;
end

end




