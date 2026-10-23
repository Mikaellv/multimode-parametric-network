%TEST_NETWORK Two small examples that check simulate_network.
%
%   Example 1: a single two-mode squeezer (2 modes), checked against the
%              analytic steady state.
%   Example 2: a two-mode squeezer (modes 1,4) feeding a beam-splitter ring
%              (modes 1,2,3) threaded by a synthetic flux Phi (4 modes). It
%              shows how to enter phases and checks gauge invariance and the
%              chirality set by the flux.
%
% Run from anywhere:  run('tests/test_network.m')   (MATLAB or Octave)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
nFail = 0;
passFail = {'PASS', 'FAIL'};
check = @(ok, msg) ok + 0*fprintf('  [%s] %s\n', passFail{2 - logical(ok)}, msg);

%% ---------------- Example 1: two-mode squeezer ----------------
fprintf('\nExample 1: two-mode squeezer between modes 1 and 2\n');
kappa = 1; g = 0.3;                                  % threshold at g = kappa/2
sim = simulate_network(2, [1 2 g], [], 'Kappa', kappa);

% Analytic: n = <a^+a> = 2g^2/(kappa^2 - 4g^2), m = <a1 a2> = kappa n/(2g)
n = 2*g^2/(kappa^2 - 4*g^2);
m = kappa*n/(2*g);
Z = diag([1 -1]);
sigmaExact = [(2*n+1)*eye(2), 2*m*Z; 2*m*Z, (2*n+1)*eye(2)];
ENexact = -log2(2*n + 1 - 2*m);
fprintf('  E_N = %.6f (analytic %.6f), purity = %.4f\n', ...
    logneg(sim.sigma, 1, 2), ENexact, purity(sim.sigma));
nFail = nFail + ~check(norm(sim.sigma - sigmaExact) < 1e-10, 'CM equals analytic steady state');
nFail = nFail + ~check(abs(logneg(sim.sigma, 1, 2) - ENexact) < 1e-10, 'E_N equals analytic value');

% The phase of the pump only rotates the correlations; E_N is unchanged.
% Same physics, three ways to give the phase:
simA = simulate_network(2, [1 2 g pi/3], []);                        % 4th column
simB = simulate_network(2, {{1, 2, g*exp(1i*pi/3)}}, []);            % complex amp
simC = simulate_network(2, [1 2 g], [], 'SqueezingPhases', pi/3);    % option
nFail = nFail + ~check(norm(simA.sigma - simB.sigma) + norm(simA.sigma - simC.sigma) < 1e-12, ...
    'phase input formats agree');
nFail = nFail + ~check(abs(logneg(simA.sigma, 1, 2) - ENexact) < 1e-10, 'E_N independent of pump phase');
nFail = nFail + ~check(norm(simA.sigma - sigmaExact) > 0.1, 'pump phase rotates the correlations');

% Heterodyne monitoring (eta = 1) purifies the state and raises E_N
simM = simulate_network(2, [1 2 g], [], 'Eta', 1);
fprintf('  monitored (eta = 1): E_N = %.4f, purity = %.6f\n', logneg(simM.sigma, 1, 2), purity(simM.sigma));
nFail = nFail + ~check(abs(purity(simM.sigma) - 1) < 1e-10, 'monitored state is pure');
nFail = nFail + ~check(logneg(simM.sigma, 1, 2) > ENexact, 'monitoring increases E_N');

% Above threshold there is no unconditional steady state
try
    simulate_network(2, [1 2 0.6], []);
    ok = false;
catch err
    ok = ~isempty(strfind(err.identifier, 'unstable'));
end
nFail = nFail + ~check(ok, 'above threshold raises an error');

%% ---------------- Example 2: squeezer + beam-splitter ring ----------------
fprintf('\nExample 2: squeezer (1,4) + beam-splitter ring 1-2-3 with flux Phi\n');
g = 0.3; h = 1/sqrt(3)*kappa;
ring = @(Phi) [1 2 h 0; 2 3 h 0; 3 1 h Phi];        % [i j amp phase]
squeezers = [1 4 g];

% A beam-splitter network alone leaves the vacuum untouched
simV = simulate_network(4, [], ring(pi/2));
nFail = nFail + ~check(norm(simV.sigma - eye(8)) < 1e-12, 'passive network keeps vacuum');

% Independent check of the drift for real couplings: H = r'*Hq*r/2 with
% squeezer g/2 (x1 p4 + p1 x4) and beam splitters h/2 (xi xj + pi pj),
% so that dr/dt = (2*Omega*Hq - kappa/2) r.
Hq = zeros(8);
Hq([1 8], [8 1]) = g/2*eye(2);  Hq([2 7], [7 2]) = g/2*eye(2);
for b = [1 2; 2 3; 3 1]'
    xi = 2*b(1) - 1; xj = 2*b(2) - 1;
    Hq([xi xj], [xj xi]) = Hq([xi xj], [xj xi]) + h/2*eye(2);
    Hq([xi xj] + 1, [xj xi] + 1) = Hq([xi xj] + 1, [xj xi] + 1) + h/2*eye(2);
end
Omega = kron(eye(4), [0 1; -1 0]);
sim0 = simulate_network(4, squeezers, ring(0));
nFail = nFail + ~check(norm(sim0.A - (2*Omega*Hq - kappa/2*eye(8))) < 1e-12, ...
    'drift matches quadrature Hamiltonian');

% Entanglement of the idler (mode 4) with each ring mode vs flux
Phis = linspace(-pi, pi, 13);
E42 = zeros(size(Phis)); E43 = E42;
for k = 1:numel(Phis)
    s = simulate_network(4, squeezers, ring(Phis(k))).sigma;
    E42(k) = logneg(s, 4, 2);
    E43(k) = logneg(s, 4, 3);
end
fprintf('  Phi/pi      : %s\n', sprintf('%6.2f ', Phis/pi));
fprintf('  E_N(4,2)    : %s\n', sprintf('%6.3f ', E42));
fprintf('  E_N(4,3)    : %s\n', sprintf('%6.3f ', E43));
nFail = nFail + ~check(max(abs(E42 - fliplr(E43))) < 1e-10, ...
    'mirror symmetry: E_N(4,2) at Phi equals E_N(4,3) at -Phi');
nFail = nFail + ~check(max(abs(E42 - E43)) > 1e-2, 'flux makes the ring chiral');

% Gauge invariance: a_j -> exp(1i*th_j) a_j shifts squeezing phases by
% th_i + th_j and beam-splitter phases by th_i - th_j. Only the loop flux
% matters, so E_N and purities must not change.
Phi = pi/2;
th = [0.3, -1.1, 2.0, 0.7];
bs = ring(Phi);
bsG = bs;  bsG(:, 4) = bs(:, 4) + th(bs(:, 1))' - th(bs(:, 2))';
sqG = [1 4 g, th(1) + th(4)];
s0 = simulate_network(4, squeezers, bs).sigma;
s1 = simulate_network(4, sqG, bsG).sigma;
pairs = nchoosek(1:4, 2);
d = 0;
for k = 1:size(pairs, 1)
    d = max(d, abs(logneg(s0, pairs(k, 1), pairs(k, 2)) - logneg(s1, pairs(k, 1), pairs(k, 2))));
    d = max(d, abs(purity(s0, pairs(k, :)) - purity(s1, pairs(k, :))));
end
nFail = nFail + ~check(d < 1e-10, 'gauge transformation leaves E_N and purities unchanged');

% Physicality of every steady state: sigma + i*Omega >= 0
nFail = nFail + ~check(min(eig(s0 + 1i*Omega)) > -1e-10, 'uncertainty principle satisfied');

% Same network, written as a cell array of tuples, monitored with eta = 0.5
simCell = simulate_network(4, {[1 4 g]}, {[1 2 h], [2 3 h], {3, 1, h, Phi}}, 'Eta', 0.5);
fprintf('  eta = 0.5: E_N(4,2) = %.4f, E_N(4,3) = %.4f, purity(all) = %.4f\n', ...
    logneg(simCell.sigma, 4, 2), logneg(simCell.sigma, 4, 3), purity(simCell.sigma));
nFail = nFail + ~check(purity(simCell.sigma) > purity(s0), 'partial monitoring increases purity');

%% ---------------- summary ----------------
if nFail == 0
    fprintf('\nAll tests passed.\n');
else
    error('test_network:failed', '%d test(s) failed.', nFail);
end

