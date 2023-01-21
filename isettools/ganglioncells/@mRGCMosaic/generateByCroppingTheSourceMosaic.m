% Method to generate the mRGCMosaic by cropping the sourceMidgetRGCMosaic
function generateByCroppingTheSourceMosaic(obj, sourceMidgetRGCMosaic)

    % Get the input cone mosaic
    obj.inputConeMosaic = sourceMidgetRGCMosaic.inputConeMosaic;
 
    if (isempty(obj.eccentricityDegs))&&(isempty(obj.sizeDegs))
        % Set the eccentricity and size
        obj.eccentricityDegs = sourceMidgetRGCMosaic.eccentricityDegs;
        obj.sizeDegs = sourceMidgetRGCMosaic.sizeDegs;

        % Get the rf positions, and cone pooling matrices
        obj.rgcRFpositionsDegs = sourceMidgetRGCMosaic.rgcRFpositionsDegs;
        obj.centerConePoolingMatrix = sourceMidgetRGCMosaic.rgcRFcenterConePoolingMatrix;
        obj.surroundConePoolingMatrix = sourceMidgetRGCMosaic.rgcRFsurroundConePoolingMatrix;
        
    else
        % Cropping must happen here
        fprintf(2,'No cropping implemented yet\n')
    end

    obj.inputConesNum = size(obj.centerConePoolingMatrix,1);
    obj.rgcsNum = size(obj.centerConePoolingMatrix,2);

end
