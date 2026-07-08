
function rfPositionsMicrons = croppedPositions(rfPositionsMicrons, eccMicrons, sizeMicrons)
    if (numel(sizeMicrons) == 1)
        sizeMicrons = sizeMicrons*[1 1];
    end
    if (numel(eccMicrons) == 1)
        eccMicrons = eccMicrons*[1 1];
    end

    xRange = eccMicrons(1) + 0.5*sizeMicrons(1)*[-1 1];
    yRange = eccMicrons(2) + 0.5*sizeMicrons(2)*[-1 1];
    
    % Crop
    xPos = squeeze(rfPositionsMicrons(:,1));
    yPos = squeeze(rfPositionsMicrons(:,2));
    idx = find(...
        (xPos >= xRange(1)) & (xPos <= xRange(2)) & ...
        (yPos >= yRange(1)) & (yPos <= yRange(2)));
    rfPositionsMicrons = rfPositionsMicrons(idx,:);
end