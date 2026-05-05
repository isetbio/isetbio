%
% 
% RGCMosaicConstructor.temporalFilterEngine.windowedOneSidedTTF(theTTF)
%
%
function theWindowedTTF = windowedOneSidedTTF(theTTF)
    
error('Are we using this?')
    theWindow = hann(2*numel(theTTF)-1);
    theOneSidedWindow = theWindow(numel(theTTF):end);
    theOneSidedWindow = reshape(theOneSidedWindow, size(theTTF));

    theWindowedTTF = theTTF .* theOneSidedWindow;
end
