function mu = purity(sigma, modes)
%PURITY Purity Tr(rho^2) = 1/sqrt(det sigma) of a Gaussian state (vacuum CM = I).
%   mu = purity(sigma)          purity of the full state
%   mu = purity(sigma, modes)   purity of the reduced state of the given modes
    if nargin > 1
        sigma = reduced_cm(sigma, modes);
    end
    mu = 1/sqrt(det(sigma));
end
