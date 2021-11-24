function spikeResp = mosaicSpikes(innerRetina)
% Convert innerRetina spikes to a matrix of number of cells by time bits.
%
% Syntax:
%   spikeResp = mosaicSpikes(innerRetina)
%
% Description:
%    Convert spikes in the innerRetina object into an NxK matrix where N is
%    the number of cells and K is the time bins.
%
% Inputs:
%    innerRetina - Object. An inner retina object.
%
% Outputs:
%    spikeResp   - Matrix. A matrix of the spike response.
%
% Optional key/value pairs:
%    None.
%

for ii = 1:4  %length(innerRetina.mosaic)
    if ii <= length(innerRetina.mosaic)
        eval(['spikesout' num2str(ii), ...
            ' = RGB2XWFormat(innerRetina.mosaic{' num2str(ii), ...
            '}.get(''spikes''));']);
    else
        eval(['spikesout' num2str(ii) ' = [];']);
    end
end

% spikesout  = RGB2XWFormat(innerRetina.mosaic{1}.get('spikes'));
% spikesout2 = RGB2XWFormat(innerRetina.mosaic{2}.get('spikes'));
% spikesout3 = RGB2XWFormat(innerRetina.mosaic{3}.get('spikes'));
% spikesout4 = RGB2XWFormat(innerRetina.mosaic{4}.get('spikes'));

timeBins = max([size(spikesout1, 2), size(spikesout2, 2), ...
    size(spikesout3, 2), size(spikesout4, 2)]);

spikesoutsm = zeros(size(spikesout1, 1) + size(spikesout2, 1) + ...
    size(spikesout3, 1) + size(spikesout4, 1), timeBins, 'uint8');
spikesoutsm(1:size(spikesout1, 1), 1:size(spikesout1, 2)) = spikesout1;
spikesoutsm(size(spikesout1, 1) + [1:size(spikesout2, 1)], ...
    1:size(spikesout2, 2)) = spikesout2;

spikesoutsm(size(spikesout1, 1) + size(spikesout2, 1) + ...
    [1:size(spikesout3, 1)], 1:size(spikesout3, 2)) = spikesout3;
spikesoutsm(size(spikesout1, 1) + size(spikesout2, 1) + ...
    size(spikesout3, 1) + [1:size(spikesout4, 1)], ...
    1:size(spikesout4, 2)) = spikesout4;

clear  spikesout1 spikesout2 spikesout3 spikesout4 

%%
spikesout = double(spikesoutsm);
pointer = 0;  % (blockNum - 1) * blocklength;
spikeResp = zeros(size(spikesoutsm, 1), ceil(size(spikesoutsm, 2) / 10));
for i = 1:size(spikesoutsm, 2) / 10
    blocksize = 10;
    endval = i * blocksize;
    if endval > size(spikesout, 2), endval = size(spikesout, 2); end
    startval = (i - 1) * blocksize + 1;
    spikeResp(:, pointer + i) = sum(spikesout(:, startval:endval), 2);
end

end
