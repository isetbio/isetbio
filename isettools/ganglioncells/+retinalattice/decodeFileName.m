function [whichEye, neuronType, sourceLatticeSizeDegs, includesProgress] = ...
        decodeFileName(theLatticeFileName)

    % Determine the eye
    if (~contains(theLatticeFileName, 'left_eye'))
        whichEye = 'right eye';
    else
        whichEye = 'left eye';
    end

    % Determine if this is a cone latice
    targetNeuronTypeString = 'cones';
    if (contains(theLatticeFileName, targetNeuronTypeString))
        neuronType = 'cones';
        neuronTypeStringFound = targetNeuronTypeString;
    end

    % Determine if this is a midget ganglion cell lattice
    targetNeuronTypeString = 'midget_ganglion_cells';
    if (contains(theLatticeFileName, targetNeuronTypeString))
        neuronType = 'midget ganglion cells';
        neuronTypeStringFound = targetNeuronTypeString;
    end

    % Determine sourceLatticeSizeDegs
    k1 = strfind(theLatticeFileName, neuronTypeStringFound);
    k2 = strfind(theLatticeFileName, 'deg');
    sourceLatticeSizeDegs = str2double(theLatticeFileName(k1+numel(neuronType)+1:k2-1));

    % Determine if this file has lattice generation progress data
    if (~contains(theLatticeFileName, 'progress'))
        includesProgress = false;
    else
        includesProgress = true;
    end

end

