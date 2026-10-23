function EN = logneg(sigma, modesA, modesB)
%LOGNEG Logarithmic negativity (base 2) between two groups of modes.
%   EN = logneg(sigma, modesA, modesB)
%   The reduced state of modesA and modesB is partially transposed on modesB
%   (p -> -p) and EN = sum_k max(0, -log2(nu_k)) over its symplectic
%   eigenvalues nu_k (vacuum CM = I). For two single modes this equals the
%   usual two-mode formula.
    s = reduced_cm(sigma, [modesA(:); modesB(:)]');
    nA = numel(modesA); nB = numel(modesB);
    P = diag([ones(1, 2*nA), repmat([1 -1], 1, nB)]);
    nu = symplectic_eigs(P*s*P);
    EN = sum(max(0, -log2(nu)));
end
