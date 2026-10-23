function s = apply_loss(sigma, transmissivity)
%APPLY_LOSS Pure-loss channel on every mode: sigma -> T*sigma + (1 - T)*I.
%   With T = cos(alpha)^2 this is the GC(P, alpha, N) channel used in
%   the original graphtestsasym.m.
    s = transmissivity*sigma + (1 - transmissivity)*eye(size(sigma));
end
