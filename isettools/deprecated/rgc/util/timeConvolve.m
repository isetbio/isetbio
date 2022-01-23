function resp = timeConvolve(mosaic, input, csFlag)
% The 1D temporal convolution of the array of center and surround responses
%
% Syntax:
%   resp = fullConvolve(mosaic, input, csFlag);
%
% Description:
%    Calculate and return the 1D temporal convolution of the array of
%    center and surround responses.
%
% Inputs: 
%    mosaic - Object. The rgc mosaic object.
%    input  - Matrix. A 4D matrix of (x, y, t, color channel)
%    csFlag - Char. A character to indicate the flag type. Options are 'c'
%             for center or 's' for surround flag.
%
% Outputs:
%    resp   - Matrix. The response over t frames.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   spConvolve.
%

% History:
%    XX/XX/15  JRG  ISETBIO TEAM, 2015
%    05/30/19  JNM  Documentation pass

%% Find bounds for size of input and output
nSamples = size(input, 3);   % Temporal samples
nChannels = size(input, 4);  % Usually just 1, but for osDisplayRGB
nCells = mosaic.get('mosaic samples'); 

%% Convolutions for each cell and channel
resp = zeros([nCells, nSamples, nChannels]);
for cc = 1:nChannels
    for ii = 1:nCells(1)
        for jj = 1:nCells(2)
            % Each cell has its own temporal impulse
            switch csFlag
                case 'c'
                    impulseR = mosaic.get('tCenter', 'cell', [ii, jj]);
                    tonicDrive = ...
                        mosaic.get('tonicDrive', 'cell', [ii, jj]);
                case 's'
                    impulseR = mosaic.get('tSurround', 'cell', [ii, jj]);
                    tonicDrive = 0;
                otherwise
            end

            % Convolve them all
            thisInput = squeeze(input(ii, jj, :, cc));
            resp(ii, jj, :, cc) = ...
                conv(thisInput, impulseR', 'same') + tonicDrive;

%             if size(thisInput, 1) > size(impulseR, 1)
%                 impulseRZP = [impulseR; ...
%                     zeros(-size(impulseR, 1) + size(thisInput, 1), 1)];
%                 thisInputZP = thisInput;
%             end
%             resp(ii, jj, :, cc) = ifft(...
%                 fft(impulseRZP) .* fft(thisInputZP)) + tonicDrive;           
        end
        % If we have multiple color channels, for the displayRGB case, sum
        % across them
        resp = sum(resp, 4);
    end
end

end