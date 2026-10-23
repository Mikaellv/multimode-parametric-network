function sim = simulate_network(N, squeezings, beamsplitters, varargin)
%SIMULATE_NETWORK Steady state of a network of N damped bosonic modes coupled
%by two-mode squeezing and beam-splitter interactions.
%
%   sim = simulate_network(N, squeezings, beamsplitters)
%   sim = simulate_network(N, squeezings, beamsplitters, 'Name', value, ...)
%
%   Couplings are given as lists of tuples. Each list may be
%     * a numeric matrix with one row per coupling:  [i j amp]  or  [i j amp phase]
%     * a cell array of tuples:  {[i j amp], [i j amp phase], {i, j, amp}, ...}
%     * [] or {} for no couplings of that kind.
%   amp may be complex; an optional phase (radians) multiplies it by exp(1i*phase).
%   Phases can alternatively be supplied as vectors with the 'SqueezingPhases'
%   and 'BeamsplitterPhases' options (one entry per tuple, added to any phase
%   already in the tuples).
%
%   Model (hbar = 1, modes a_1..a_N, rotating frame):
%     squeezing (i,j,g):     H += 1i*(g a_i^+ a_j^+ - g^* a_i a_j)  ->  da_i/dt += g a_j^+
%                            (i == j gives single-mode squeezing, da_i/dt += g a_i^+)
%     beam splitter (i,j,h): H += h a_i^+ a_j + h^* a_j^+ a_i       ->  da_i/dt += -1i h a_j
%     every mode is damped at rate kappa_j into a vacuum bath.
%   Quadratures x = a + a^+, p = -1i(a - a^+), ordered (x1,p1,x2,p2,...), so
%   the vacuum covariance matrix (CM) is the identity.
%
%   Options:
%     'Kappa'              damping rate(s), scalar or 1xN (default 1)
%     'Eta'                heterodyne efficiency of the monitored output of each
%                          mode, scalar or 1xN in [0,1] (default 0). Eta = 0 gives
%                          the unconditional steady state (Lyapunov equation);
%                          Eta > 0 gives the conditional steady state under
%                          continuous heterodyne monitoring (Riccati equation).
%     'Detunings'          mode frequencies in the rotating frame, 1xN (default 0)
%     'SqueezingPhases'    extra phase for each squeezing tuple (default 0)
%     'BeamsplitterPhases' extra phase for each beam-splitter tuple (default 0)
%
%   Output struct fields:
%     sigma   2N x 2N steady-state covariance matrix (vacuum = identity)
%     A, D    quadrature drift and diffusion matrices (dsigma/dt = A s + s A' + D - ...)
%     M       complex drift in the (a; a^+) basis
%     S, Hbs  symmetric squeezing matrix and Hermitian beam-splitter matrix
%     stable  true if the unconditional dynamics is below threshold
%     residual  norm of the residual of the solved steady-state equation
%
%   See also LOGNEG, PURITY, REDUCED_CM, APPLY_LOSS.

    p = inputParser;
    p.addParameter('Kappa', 1);
    p.addParameter('Eta', 0);
    p.addParameter('Detunings', zeros(1, N));
    p.addParameter('SqueezingPhases', []);
    p.addParameter('BeamsplitterPhases', []);
    p.parse(varargin{:});
    opt = p.Results;

    kappa = expand(opt.Kappa, N, 'Kappa');
    eta = expand(opt.Eta, N, 'Eta');
    delta = expand(opt.Detunings, N, 'Detunings');
    if any(eta < 0 | eta > 1)
        error('simulate_network:eta', 'Eta must lie in [0, 1].');
    end

    sq = parse_couplings(squeezings, opt.SqueezingPhases, N, 'squeezing');
    bs = parse_couplings(beamsplitters, opt.BeamsplitterPhases, N, 'beam splitter');

    % Symmetric squeezing matrix S and Hermitian beam-splitter matrix Hbs
    S = zeros(N);
    for k = 1:size(sq, 1)
        i = sq(k, 1); j = sq(k, 2); g = sq(k, 3);
        S(i, j) = S(i, j) + g;
        if i ~= j
            S(j, i) = S(j, i) + g;
        end
    end
    Hbs = diag(delta);
    for k = 1:size(bs, 1)
        i = bs(k, 1); j = bs(k, 2); h = bs(k, 3);
        if i == j
            error('simulate_network:bs', ...
                'Beam splitter (%d,%d) couples a mode to itself; use Detunings.', i, j);
        end
        Hbs(i, j) = Hbs(i, j) + h;
        Hbs(j, i) = Hbs(j, i) + conj(h);
    end

    % Heisenberg-Langevin drift in the (a; a^+) basis, transformed to quadratures
    Kap = diag(kappa);
    M = [-1i*Hbs - Kap/2, S; conj(S), 1i*conj(Hbs) - Kap/2];
    T = quadrature_transform(N);
    A = real(T*M*T'/2);                  % T^-1 = T'/2
    D = kron(Kap, eye(2));

    stable = all(real(eig(A)) < 0);

    if all(eta == 0)
        if ~stable
            error('simulate_network:unstable', ...
                'Network is above threshold: no unconditional steady state.');
        end
        sigma = sylvester(A, A', -D);
        sigma = (sigma + sigma')/2;
        residual = norm(A*sigma + sigma*A' + D, 'fro');
    else
        % Heterodyne monitoring: A s + s A' + D - (s - I) W (s - I) = 0
        W = kron(diag(eta.*kappa/2), eye(2));
        At = A + W;
        sigma = solve_riccati(At, W, D - W);
        I = eye(2*N);
        residual = norm(A*sigma + sigma*A' + D - (sigma - I)*W*(sigma - I), 'fro');
    end
    if residual > 1e-8*max(1, norm(sigma, 'fro'))
        warning('simulate_network:residual', ...
            'Steady-state residual is large (%.2e).', residual);
    end

    sim = struct('sigma', sigma, 'A', A, 'D', D, 'M', M, 'S', S, 'Hbs', Hbs, ...
        'kappa', kappa, 'eta', eta, 'N', N, 'stable', stable, 'residual', residual);
end

function v = expand(v, N, name)
    if isscalar(v)
        v = v*ones(1, N);
    elseif numel(v) == N
        v = reshape(v, 1, N);
    else
        error('simulate_network:size', '%s must be a scalar or have N entries.', name);
    end
end

function c = parse_couplings(list, extraPhases, N, kind)
    % Returns rows [i j complexAmplitude]
    if isempty(list)
        c = zeros(0, 3);
        rows = {};
    elseif iscell(list)
        rows = cellfun(@(t) cell2mat_safe(t), list, 'UniformOutput', false);
    else
        rows = num2cell(list, 2);
    end
    n = numel(rows);
    if isempty(extraPhases)
        extraPhases = zeros(1, n);
    elseif numel(extraPhases) ~= n
        error('simulate_network:phases', ...
            'Number of %s phases (%d) does not match number of tuples (%d).', ...
            kind, numel(extraPhases), n);
    end
    if n > 0
        c = zeros(n, 3);
    end
    for k = 1:n
        t = rows{k};
        if numel(t) < 3 || numel(t) > 4
            error('simulate_network:tuple', ...
                'Each %s tuple must be (i, j, amp) or (i, j, amp, phase).', kind);
        end
        i = t(1); j = t(2);
        if any([i j] ~= round(real([i j]))) || any([i j] < 1) || any([i j] > N)
            error('simulate_network:index', ...
                '%s tuple %d has mode indices outside 1..%d.', kind, k, N);
        end
        ph = extraPhases(k);
        if numel(t) == 4
            ph = ph + real(t(4));
        end
        c(k, :) = [real(i), real(j), t(3)*exp(1i*ph)];
    end
end

function v = cell2mat_safe(t)
    if iscell(t)
        v = [t{:}];
    else
        v = t;
    end
end

function T = quadrature_transform(N)
    % (x1,p1,...,xN,pN)^T = T * (a_1..a_N, a_1^+..a_N^+)^T
    T = zeros(2*N);
    for j = 1:N
        T(2*j-1, j) = 1;   T(2*j-1, N+j) = 1;
        T(2*j,   j) = -1i; T(2*j,   N+j) = 1i;
    end
end

function X = solve_riccati(F, W, Q)
    % Stabilizing solution of F X + X F' - X W X + Q = 0 via the ordered
    % Schur form of the Hamiltonian matrix (no toolbox needed).
    n = size(F, 1);
    Z = [F', -W; -Q, -F];
    [U, R] = schur(Z, 'complex');
    [U, ~] = ordschur(U, R, real(diag(R)) < 0);
    X = U(n+1:end, 1:n)/U(1:n, 1:n);
    X = real(X + X')/2;
end
