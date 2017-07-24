function resp = timeConvolve(mosaic, input, csFlag)
% The 1D temporal convolution of the array of center and surround responses
%
%   resp = fullConvolve(mosaic, input, csFlag);
%
% Inputs: 
%  rgc mosaic   - 
%  respCenter   - 
%  respSurround -
%   4D matrix of (x,y,t,color channel)
%  csFlag - Center or surround flag ('c' or 's')
%
% Outputs: 
%   The response over t frames
%
% See also: spConvolve.
%
% JRG, ISETBIO TEAM, 2015

%% Find bounds for size of input and output

nSamples    = size(input,3);    % Temporal samples
nChannels = size(input,4);    % Usually just 1, but for osDisplayRGB
nCells      = mosaic.get('mosaic samples'); 

%% Convolutions for each cell and channel

resp = zeros([nCells,nSamples,nChannels]);

for cc = 1:nChannels
    for ii = 1:nCells(1)
        for jj = 1:nCells(2)
            
            % Each cell has its own temporal impulse
            switch csFlag
                case 'c'
                    impulseR = mosaic.get('tCenter','cell',[ii,jj]);
                    tonicDrive = mosaic.get('tonicDrive','cell',[ii,jj]);
                case 's'
                    impulseR = mosaic.get('tSurround','cell',[ii,jj]);
                    tonicDrive = 0;
                otherwise
            end
            
            % Convolve them all
            thisInput = squeeze(input(ii,jj,:,cc));
            resp(ii,jj,:,cc) = conv(thisInput, impulseR','same')+tonicDrive;
            
%             if size(thisInput,1) > size(impulseR,1)
%                 impulseRZP = [impulseR; zeros(-size(impulseR,1)+size(thisInput,1),1)];
%                 thisInputZP = thisInput;
%             end
%             resp(ii,jj,:,cc) = ifft(fft(impulseRZP).*fft(thisInputZP))+tonicDrive;           
            
        end
        
        % If we have multiple color channels, for the displayRGB case, sum
        % across them
        resp = sum(resp, 4);
    end
end

end