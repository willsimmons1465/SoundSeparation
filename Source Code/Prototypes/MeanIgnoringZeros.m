function [ result ] = MeanIgnoringZeros( A )
%MeanIgnoringZeros Takes the mean of the non-zero elements of a row

S = sum(A')';
N = sum(A' ~= 0)';
result = S./N;

end

