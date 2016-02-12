function [SpRast, binctrs] = binSpTimes(sptimes, bins, slen);
%  [Sp, binctrs] = binSpTimes(sptimes, bins, slen);
%  Convert floating point spike times into a binned representation
%
%  Inputs:  sptimes = vector or matrix of spike times
%             If a matrix, each column is the spike times for a repeat 
%             (padded with zeros)
%           bins =  either the bin width, or a vector of bin centers
%           slen = length of desired spike raster (OPTIONAL)
%
%  Outputs:  Sp = binned spike train
%            binctrs = centers of spike train time bins 


% Check that sptimes are monotonic increasing
% Removed: fails on matrices
%isi = diff(sptimes);
%if min(isi) <= 0
%    fprintf(1, ['WARNING: spike times NOT monotonically increasing!! ' ...
%                '(binSpTimes)\n']);
%    keyboard;
%end

if length(bins) == 1
    if nargin == 2
        slen = ceil(max(max(sptimes))/bins)*bins;
    end
    binctrs = [-bins/2:bins:(slen+bins)]';
    SpRast = hist(sptimes, binctrs);  
    
    if size(SpRast,1) == 1  %  If 1D, make into column vector;
        SpRast = SpRast';
    end

    SpRast = SpRast(2:end-1,:);
    binctrs = binctrs(2:end-1);
    
else

    binctrs = bins;
    bins = bins(2)-bins(1);
    SpRast = hist(sptimes, binctrs);
    
end

