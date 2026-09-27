function run_all
% Runs the pipeline of "Electromagnetic field imaging in arbitrary scattering
% environments" (K. Sastry, C. Bhat, R. Solimene, U. K. Khankhoje, IEEE TCI 2021)
% for the paper's 4-object room at a sampling rate of 0.3 (194 measurements) and
% 25 dB SNR, once with MPME sensor locations and once with random locations.
% Takes about 6-7 minutes on a laptop.
%
% Requires CVX (http://cvxr.com/cvx) on the MATLAB path.
% Intermediate .mat files are written to ./results, figures to ./figures.

root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'code'));
out = fullfile(root, 'results');
if ~exist(out, 'dir'), mkdir(out); end
here = pwd; cleanup = onCleanup(@() cd(here));
cd(out);

% Run 1: geometry, forward solve, true fields, MPME locations, CS-SOM, prediction matrix
cfg = struct('run_obj',true,'run_fwd',true,'run_grd',true,'gen_invp_spec',true, ...
    'run_pts',false,'run_pts_mpme',true,'run_inv',true,'run_ttf',true,'run_prd',true, ...
    'plot_2d',false,'scheme','mpme');
fprintf('\n=== Run 1/2: MPME sensor placement\n');
run_fmap(cfg);

% Run 2: same problem with random locations (reuses the matrices from run 1)
cfg = struct('run_obj',false,'run_fwd',false,'run_grd',false,'gen_invp_spec',false, ...
    'run_pts',true,'run_pts_mpme',false,'run_inv',true,'run_ttf',false,'run_prd',false, ...
    'plot_2d',false,'scheme','random');
fprintf('\n=== Run 2/2: random sensor placement\n');
run_fmap(cfg);

make_figures(out, fullfile(root, 'figures'));
end

function run_fmap(cfg) %#ok<INUSD>
% fmap_data_3 reads its stage flags from cfg (it clears everything else)
fmap_data_3
close all
end
