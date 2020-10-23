function parseInput(neuronType, whichEye)
    validNeuronTypes = {'cones', 'midget ganglion cells'};
    validEyes = {'left eye', 'right eye'};
    assert(ismember(neuronType, validNeuronTypes), sprintf('Unknown neuron type: ''%s''.', neuronType));
    assert(ismember(whichEye, validEyes), sprintf('Unknown eye: ''%s''.', whichEye));
end

