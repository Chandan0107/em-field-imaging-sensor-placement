function make_figures(results_dir, fig_dir)
% Figures from the two runs of fmap_data_3.m: sensor locations (cf. paper Fig. 3),
% true vs predicted total field (Fig. 7) and relative error, MPME vs random (Fig. 10).

if nargin < 1, results_dir = 'results'; end
if nargin < 2, fig_dir = 'figures'; end
if ~isfolder(fig_dir), mkdir(fig_dir); end
M = load(fullfile(results_dir, 'result_mpme.mat'));
R = load(fullfile(results_dir, 'result_random.mat'));
lambda = M.invp.lambda;
bl = M.invp.corners(1).bl; tr = M.invp.corners(1).tr; c = (bl + tr) / 2;

toimg = @(S, v) to_image(S.locs_gp, v, c, lambda);
true_field = M.sol_grid(M.ind) + M.sol_gridi(M.ind);
[Et, ax, ay] = toimg(M, true_field);
Em = toimg(M, M.est_grid);
Er = toimg(R, R.est_grid);
relerr = @(E) min(abs(E - Et) ./ abs(Et), 1);

%% Sensor locations
fig1 = figure('Color', 'w', 'Position', [100 100 1000 460]);
S = {M, R}; names = {'MPME', 'Random'};
for p = 1:2
    subplot(1, 2, p); hold on; box on
    draw_objects(M.invp, c, lambda);
    xy = (S{p}.locs - c) / lambda;
    plot(xy(:, 1), xy(:, 2), 'k.', 'MarkerSize', 10);
    axis equal; axis([ax(1) ax(end) ay(1) ay(end)] + [-0.3 0.3 -0.3 0.3]);
    xlabel('x [\lambda]'); ylabel('y [\lambda]');
    title(sprintf('%s: %d measurements', names{p}, size(S{p}.locs, 1))); set(gca, 'FontSize', 12);
end
exportgraphics(fig1, fullfile(fig_dir, 'sampling_locations.png'), 'Resolution', 150);

%% True vs predicted field, and relative errors
fig2 = figure('Color', 'w', 'Position', [100 100 1400 380]);
cl = [0 prctile(abs(Et(:)), 99)];
panels = {abs(Et), abs(Em), relerr(Em), relerr(Er)};
titles = {'True total field |E|', sprintf('Predicted, MPME (\\Delta_G = %.1f%%)', 100*M.gerror), ...
          sprintf('Relative error, MPME (%.1f%%)', 100*M.gerror), ...
          sprintf('Relative error, random (%.1f%%)', 100*R.gerror)};
for p = 1:4
    subplot(1, 4, p);
    h = imagesc(ax, ay, panels{p}); set(h, 'AlphaData', ~isnan(panels{p}));
    axis xy image; colorbar; set(gca, 'Color', [0.85 0.85 0.85]);
    if p <= 2, caxis(cl); colormap(gca, turbo); else, caxis([0 1]); colormap(gca, hot); end
    xlabel('x [\lambda]'); ylabel('y [\lambda]'); title(titles{p}); set(gca, 'FontSize', 11);
end
exportgraphics(fig2, fullfile(fig_dir, 'field_prediction.png'), 'Resolution', 150);
fprintf('Figures saved to %s\n', fig_dir);
end

function [img, ax, ay] = to_image(locs, v, c, lambda)
% Place values given at grid points back onto a regular image (NaN where absent).
xs = unique(round(locs(:, 1), 10)); ys = unique(round(locs(:, 2), 10));
[~, ix] = ismember(round(locs(:, 1), 10), xs);
[~, iy] = ismember(round(locs(:, 2), 10), ys);
img = nan(numel(ys), numel(xs));
img(sub2ind(size(img), iy, ix)) = v;
ax = (xs - c(1)) / lambda; ay = (ys - c(2)) / lambda;
end

function draw_objects(invp, c, lambda)
% Approximate object contours used by the solver (the wall is surface 1).
for s = 2:numel(invp.surf)
    m = (invp.surf(s).m - c) / lambda;
    fill(m(:, 1), m(:, 2), [0.95 0.55 0.5], 'EdgeColor', 'none');
end
w = (invp.surf(1).m - c) / lambda;
plot(w([1:end 1], 1), w([1:end 1], 2), 'b-', 'LineWidth', 2);
end
