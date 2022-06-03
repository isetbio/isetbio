function result = current(obj,excitations,irf,timeAxis)
% Calculate the photocurrent from the IRFs and excitations
%
% Inputs
%  excitations - nInstances, nTimes, nCones
%  irf - Three impulse response functions
%  timeAxis
%
% Output
%  photocurrent time series
%
% See also
%  currentIRF
%  t_cMosaicLetter for an example usage

result = zeros(size(excitations));
nInstances = size(excitations,1);

for ii=1:nInstances
    idx = obj.lConeIndices;
    tmp = squeeze(excitations(ii,:,idx));
    for jj=1:numel(idx)
        result(ii,:,idx(jj)) = cconv(tmp(:,jj),irf(:,1),numel(timeAxis));
    end   
end
for ii=1:nInstances
    idx = obj.mConeIndices;
    tmp = squeeze(excitations(ii,:,idx));
    for jj=1:numel(idx)
        result(ii,:,idx(jj)) = cconv(tmp(:,jj),irf(:,2),numel(timeAxis));
    end   
end
for ii=1:nInstances
    idx = obj.sConeIndices;
    tmp = squeeze(excitations(ii,:,idx));
    for jj=1:numel(idx)
        result(ii,:,idx(jj)) = cconv(tmp(:,jj),irf(:,3),numel(timeAxis));
    end   
end

end