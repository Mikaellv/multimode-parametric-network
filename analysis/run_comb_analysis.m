%RUN_COMB_ANALYSIS Reproduces the analysis of graphtests.m / graphtestsasym.m
%with simulate_network.
%
% A frequency comb of 44 modes (labels +-1..+-22) is pumped by 1..15 pumps of
% equal amplitude r. For each number of pumps we compute the conditional
% steady state under continuous heterodyne monitoring of all outputs and
% evaluate the logarithmic negativity E_N and the purity of modes (-1, 1).
%
%   'sym'  : pump k couples comb pairs with frequency sum s_k = 0, 2, -2, 4, ...
%   'asym' : pumps 1-2 as in 'sym', every further pump couples each comb mode
%            to its own idler mode outside the comb (44 new modes per pump)
%
% Six cases as in the original figure: {sym, asym} x {equal phases, random
% phases, random phases + noise}. "Noise" is heterodyne efficiency eta = 1/2,
% i.e. b = pi/4 in the original sec(b)^2 + n*tan(b)^2 heterodyne noise with n = 1.
% Values hardcoded in graphtests.m are overlaid as markers for comparison;
% random-phase curves only agree qualitatively there, since the original
% phase draws were unseeded.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---------------- parameters ----------------
r = 0.25;             % pump amplitude (units of the damping rate kappa = 1)
nPairs = 22;          % comb modes -22..-1, 1..22
nPumpsMax = 15;
etaIdeal = 1;         % ideal heterodyne
etaNoisy = 0.5;       % "noise" case
seed = 1;             % random-phase seed
modesOfInterest = comb_index([-1, 1]);

rng(seed);
thetaRand = pi*(1 - 2*rand(1, nPumpsMax));
amps = r*ones(1, nPumpsMax);

cases = struct( ...
    'geometry', {'sym', 'sym', 'sym', 'asym', 'asym', 'asym'}, ...
    'phases',   {zeros(1, nPumpsMax), thetaRand, thetaRand, ...
                 zeros(1, nPumpsMax), thetaRand, thetaRand}, ...
    'eta',      {etaIdeal, etaIdeal, etaNoisy, etaIdeal, etaIdeal, etaNoisy}, ...
    'label',    {'Sym equal phases no noise', 'Sym rand phases no noise', ...
                 'Sym rand phases noise', 'Asym equal phases no noise', ...
                 'Asym rand phases no noise', 'Asym rand phases noise'});

% ---------------- simulation ----------------
EN = zeros(numel(cases), nPumpsMax);
PUR = zeros(numel(cases), nPumpsMax);
for c = 1:numel(cases)
    for np = 1:nPumpsMax
        nModes = 2*nPairs;
        if strcmp(cases(c).geometry, 'asym') && np >= 3
            nModes = 2*nPairs*(np - 1);
        end
        sq = comb_pumps(np, nPairs, amps, cases(c).phases, cases(c).geometry);
        sim = simulate_network(nModes, sq, [], 'Eta', cases(c).eta);
        EN(c, np) = logneg(sim.sigma, modesOfInterest(1), modesOfInterest(2));
        PUR(c, np) = purity(sim.sigma, modesOfInterest);
    end
    fprintf('%-28s E_N   = %s\n', cases(c).label, sprintf('%.4f ', EN(c, :)));
    fprintf('%-28s Purity= %s\n', '', sprintf('%.4f ', PUR(c, :)));
end

% Values hardcoded in graphtests.m (lognegs1..3, lognegsas1..3, purs..., pursas...)
origEN = [ ...
    0.6942, 0.4820, 0.3369, 0.2334, 0.1556, 0.0944, 0.0447, 0.0032, 0, 0, 0, 0, 0, 0, 0;
    0.6942, 0.4820, 0.3030, 0.1754, 0.0615, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0.6309, 0.3874, 0.1653, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    NaN, 0.4820, 0.3061, 0.1571, 0.0284, NaN, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0.6942, 0.4820, 0.3061, NaN, 0.0284, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0.6309, 0.3874, 0.1702, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
origPUR = [ ...
    1.0000, 0.8139, 0.7049, 0.6352, 0.5869, 0.5513, 0.5237, 0.5017, 0.4839, 0.4686, 0.4526, 0.4380, 0.4270, 0.4185, 0.4113;
    1.0000, 0.8139, 0.6814, 0.5989, 0.5329, 0.4631, 0.4267, 0.4034, 0.3870, 0.3554, 0.3228, 0.3142, 0.2984, 0.2732, 0.2653;
    0.8941, 0.6616, 0.5137, 0.3972, 0.3326, 0.2908, 0.2531, 0.2228, 0.1946, 0.1699, 0.1552, 0.1445, 0.1320, 0.1209, 0.1118;
    NaN, NaN, 0.6835, 0.5878, 0.5149, 0.4577, 0.4117, 0.3740, NaN, 0.3158, 0.2930, 0.2732, NaN, 0.2406, 0.2270;
    1.0000, 0.8139, 0.6835, 0.5878, 0.5149, 0.4577, NaN, 0.3740, 0.3425, 0.3158, NaN, 0.2732, 0.2558, 0.2406, NaN;
    0.8941, 0.6616, 0.5114, 0.4107, 0.3359, 0.2821, 0.2413, 0.2097, 0.1846, 0.1644, 0.1478, 0.1339, 0.1223, 0.1123, 0.1038];
for c = [1 4]   % deterministic cases: compare directly
    dEN = abs(EN(c, :) - origEN(c, :));
    dPUR = abs(PUR(c, :) - origPUR(c, :));
    fprintf('max |E_N - original| (%s): %.1e,  max |purity - original|: %.1e\n', ...
        cases(c).label, max(dEN(~isnan(dEN))), max(dPUR(~isnan(dPUR))));
end

% ---------------- save ----------------
pumps = 1:nPumpsMax;
save(fullfile(outdir, 'comb_results.mat'), 'pumps', 'EN', 'PUR', 'cases', ...
    'thetaRand', 'r', 'etaIdeal', 'etaNoisy', 'seed', '-v7');
fid = fopen(fullfile(outdir, 'comb_results.csv'), 'w');
fprintf(fid, 'case,quantity%s\n', sprintf(',pumps_%d', pumps));
for c = 1:numel(cases)
    fprintf(fid, '%s,E_N%s\n', cases(c).label, sprintf(',%.6f', EN(c, :)));
    fprintf(fid, '%s,purity%s\n', cases(c).label, sprintf(',%.6f', PUR(c, :)));
end
fclose(fid);

% ---------------- plots (style of graphtests.m) ----------------
colors = {[1 0 0], [0 0.5 0], [1 0 1], [0 0 1], [0.7 0.6 0], [0 0.8 0.8]};
styles = {'-', '-', '-', ':', ':', ':'};
lw = 3;
quantities = {EN, PUR; origEN, origPUR};
ylabels = {'E_N', 'Purity'};
files = {'symasymcomplogneg.png', 'symasymcomppur.png'};
for q = 1:2
    figure(q); clf; hold on
    h = [];
    for c = 1:numel(cases)
        h = [h, plot(pumps, quantities{1, q}(c, :), 'Color', colors{c}, ...
            'LineStyle', styles{c}, 'LineWidth', lw)]; %#ok<AGROW>
        plot(pumps, quantities{2, q}(c, :), 'o', 'Color', colors{c}, 'MarkerSize', 6);
    end
    hold off
    xlim([1 nPumpsMax]);
    if q == 2, ylim([0 1]); end
    set(gca, 'XTick', 1:2:nPumpsMax, 'FontSize', 14);
    xlabel('Number of pumps', 'FontSize', 16);
    ylabel(ylabels{q}, 'FontSize', 16);
    legend(h, {cases.label}, 'Location', 'northeast');
    title('lines: simulate\_network,  circles: graphtests.m', 'FontSize', 12);
    grid on
    print(gcf, fullfile(outdir, files{q}), '-dpng', '-r150');
end

