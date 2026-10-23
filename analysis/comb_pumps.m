function sq = comb_pumps(nPumps, nPairs, amps, phases, geometry)
%COMB_PUMPS Squeezing tuples [i j amp phase] for the pumped frequency comb of
%graphtests.m ('sym') and graphtestsasym.m ('asym').
%
%   The comb has 2*nPairs modes with labels l = +-1..+-nPairs at frequency
%   offsets f(l) = l - sign(l)/2 around half the pump frequency.
%   'sym':  pump k couples every pair of comb modes with f(l1)+f(l2) = s_k,
%           s_k = 0, 2, -2, 4, -4, ... for k = 1, 2, 3, ...
%   'asym': pumps 1 and 2 as above; every further pump k couples each comb
%           mode one-to-one to its own fresh idler mode (indices
%           2*nPairs*(k-2) + 1 .. 2*nPairs*(k-1)).
%   amps, phases: 1 x nPumps amplitude and phase of each pump.
    if nargin < 5, geometry = 'sym'; end
    labels = [-nPairs:-1, 1:nPairs];
    f = labels - sign(labels)/2;
    nComb = 2*nPairs;
    sq = zeros(0, 4);
    for k = 1:nPumps
        if strcmp(geometry, 'asym') && k >= 3
            m = (1:nComb)';
            sq = [sq; m, nComb*(k-2) + m, repmat([amps(k), phases(k)], nComb, 1)]; %#ok<AGROW>
            continue
        end
        if k == 1, s = 0; elseif mod(k, 2) == 0, s = k; else, s = -(k-1); end
        for a = 1:numel(labels)
            for b = a+1:numel(labels)
                if f(a) + f(b) == s
                    sq(end+1, :) = [comb_index(labels(a)), comb_index(labels(b)), ...
                        amps(k), phases(k)]; %#ok<AGROW>
                end
            end
        end
    end
end
