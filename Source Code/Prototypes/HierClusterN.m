function [ Mapping ] = HierClusterN( D, N )
%HierClusterN - Split into N clusters by Hierarchical Clustering
%   Inputs: D is the distance matrix of samples, N is the number of desired
%   clusters
%   Outputs: Mapping is a vector with each index being a value from 0 to
%   N-1 indicating the cluster in which the corresponding sample resides

%Handle Base Case
NumOfSamples = int32(size(D,1));
if NumOfSamples <= N
    Mapping = 0:NumOfSamples-1;
    return;
end

%Ignore Diagonals
for i = 1:NumOfSamples
    D(i,i) = inf;
end

%Identify candidates for merging (should give lower index as
%PairToMerge(1))
MinLocs = int32(find(D == min(min(D))));
PairToMerge = [idivide(MinLocs(1)-1, NumOfSamples, 'floor')+1, mod(MinLocs(1)-1, NumOfSamples)+1];

if PairToMerge(1) > PairToMerge(2)
    PairToMerge = fliplr(PairToMerge);
end

%Perform merge on matrix
DMod = D;
DMod(PairToMerge(1),:) = 0.5*(D(PairToMerge(1),:) + D(PairToMerge(2),:));
DMod(:,PairToMerge(1)) = 0.5*(D(:,PairToMerge(1)) + D(:,PairToMerge(2)));
DReduced = DMod([1:PairToMerge(2)-1, PairToMerge(2)+1:end],[1:PairToMerge(2)-1, PairToMerge(2)+1:end]);

RecMapping = HierClusterN(DReduced, N);
if PairToMerge(2) == NumOfSamples
    Mapping = [RecMapping(1:end) RecMapping(PairToMerge(1))];
else
    Mapping = [RecMapping(1:PairToMerge(2)-1) RecMapping(PairToMerge(1)) RecMapping(PairToMerge(2):end)];
end

end

