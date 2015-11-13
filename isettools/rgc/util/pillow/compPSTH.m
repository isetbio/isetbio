function [psth, tt, pstv, spRaster] = compPSTH(Mtsp, dt, sigma, trnge, sig2, type);
%  [psth, tt, pstv, spRaster] = compPSTH(Mtsp, dt, sigma, trnge, sig2, type);
%  
%  Compute PSTH and PSTV from repeated responses
%
%  Inputs:  Mtsp = spike times (each col = one repeat).
%           dt = size of bins for computing raster
%           sigma = stdev of Gaussian for filtering PSTH
%           trnge = [strt stp];
%           sig2 = optional: sigma for variance curves
%           type = 's' = square bins,  'g' = gaussian bins
%
%  Outputs:  psth:  firing rate vs. time, avged over trials.
%            pstv:  variance over trials (square or gaussian bins)
%            pstv2:  variance over trials (square bins)
%            tt:  times corresponding to psth
%
% Sample call:  
% [psth,tt,pstv,spr] = compPSTH(Mtsp*dt, .001, .002, [0 1], .005)

if iscell(Mtsp)
    Mtsp = convrtTspCelltoMatrix(Mtsp);
end

nrpts = size(Mtsp,2);
if nargin > 3
    slen = trnge(2);
else
    slen = ceil(max(max(Mtsp)));
end

% Convert repeat spike times into rasters
[spRaster, tt] = binSpTimes(Mtsp, dt, slen); 

if nargin > 2
    ii = find((tt<=trnge(2)) & (tt>= trnge(1)));
    spRaster = spRaster(ii,:)';
    tt = tt(ii);
else
    spRaster = spRaster';
end

% Compute PSTH]
psthraw = sum(spRaster, 1)'/nrpts/dt;

% Filter w/ Gaussian
xfilt = dt*(-ceil(4*sigma/dt):1:ceil(4*sigma/dt))';
ff = ieNormpdf(xfilt, 0, sigma)*dt;

psth = conv2(psthraw, ff, 'same');

if nargout > 2
    
    if (nargin >= 5)
        sigma = sig2;
        xfilt = dt*(-ceil(4*sigma/dt):1:ceil(4*sigma/dt))';
    end
                        
    if nargin < 6
%         fprintf(1, 'compPSTH: using square bins for PSTV\n');
        type = 's';
    end
    
    if type(1) == 's'
        ff = (xfilt>=-sigma/2) & (xfilt<=sigma/2);
    elseif type(1) == 'g'
        ff = ieNormpdf(xfilt, 0, sigma)*dt;
    else
        error('compPSTH:  unknown type of bins for PSTV:  %s', type);
    end

    smrast = conv2(spRaster, double(ff'), 'same');
    pstv = var(smrast);
        
end



% 
% 
% if 1  %  Make figure;
%     if nargin < 3;
%         nsecs = 1;  %  Number of seconds to show;
%     end
%     T = floor(nsecs/dt2);  % length of vector corresponding to nsecs of stim
%     tt = (1:slen*scale)*dt2-dt2/2;
% 
%     subplax([3 1], [2 1], [2 1]);
%     imagesc([0 nsecs], [1 nrpts], -spRaster(:,1:T));
%     ylabel('Repeat #'); colormap gray;
%     
%     subplot(313);
% 
%     plot(tt(1:T), psth(1:T));
% %    hold on; plot(tt2, smoothPSTH(1:T2), 'g'); hold off;
%     ylabel('firing rate (sp/s)');
%     xlabel('time (sec)');
% end
