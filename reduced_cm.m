function s = reduced_cm(sigma, modes)
%REDUCED_CM Covariance matrix of a subset of modes (quadrature order x1,p1,...).
%   s = reduced_cm(sigma, modes)
    idx = reshape([2*modes(:)' - 1; 2*modes(:)'], 1, []);
    s = sigma(idx, idx);
end
