function nu = symplectic_eigs(sigma)
%SYMPLECTIC_EIGS Symplectic eigenvalues of a 2n x 2n covariance matrix (sorted).
    n = size(sigma, 1)/2;
    Omega = kron(eye(n), [0 1; -1 0]);
    ev = sort(abs(eig(1i*Omega*sigma)));
    nu = ev(1:2:end);       % eigenvalues come in +-nu pairs
end
