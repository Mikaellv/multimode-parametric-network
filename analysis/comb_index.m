function idx = comb_index(labels)
%COMB_INDEX Mode index of frequency-comb labels -1,1,-2,2,... (interleaved order,
%the node order MATLAB's graph() produced in the original graphtests.m).
    idx = 2*abs(labels) - (labels < 0);
end
