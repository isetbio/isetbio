%
%PAL_AMRF_pdfDescriptives   Calculate some summary statistics for a 
%   (discretized) probability density function.
%   
%   syntax: [Mode Mean SD] = PAL_AMRF_pdfDescriptives(pdf, values) 
%
%   returns the mode, mean, and standard deviation based on the densities 
%       defined in vector 'pdf' across the values defined in vector 
%       'values' (which should be equal in size to vector 'pdf').
%
%   Example: 
%
%   [Md M SD] = PAL_AMRF_pdfDescriptives([.2 .5 .25 .05], [3 4 5 6])
%   
%   returns:
%
%   Md = 4
%
%   M =  4.1500
%   
%   SD = 0.7921
%       
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [Mode, Mean, SD] = PAL_AMRF_pdfDescriptives(pdf, values)

pdf = pdf./sum(pdf);
if min(pdf) == max(pdf)                     %if uniform distribution
    Mode = values(round(length(values)/2)); %define mode as median
                                            %(assumes values are ordered)
else
    [trash, I] = max(pdf);
    Mode = values(I);
end
Mean = sum(pdf.*values);
SD = (sum(pdf.*(values-Mean).^2))^.5;