function test_ConeToMRGCMosaicConnector()

    sourceLatticeSizeDegs = 60;
    customDegsToMMsConversionFunction = @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x);
    customMMsToDegsConversionFunction = @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x);


    % Generate the input cone mosaic (the source)
    eccDegs = [0 0]; sizeDegs = 0.5*[0.5 0.35];
    %eccDegs = [5 4]; sizeDegs = 1.5*[0.5 0.35];

    theInputConeMosaic = cMosaic(...
       'sourceLatticeSizeDegs', sourceLatticeSizeDegs, ...
       'eccentricityDegs', eccDegs, ...
       'sizeDegs', sizeDegs, ...
       'coneDensities', [0.6 0.3 0.1], ...
       'overlappingConeFractionForElimination', 0.5, ...
       'rodIntrusionAdjustedConeAperture', true, ...
       'coneApertureModifiers', struct('smoothLocalVariations', true));

    
    % Source lattice (i.e. cone mosaic lattice) meta data, here cone aperture size and cone types
    metaData.coneApertureDiametersMicrons = theInputConeMosaic.coneApertureDiametersMicrons;    
    metaData.coneTypes = theInputConeMosaic.coneTypes;
    metaData.coneTypeIDs = [theInputConeMosaic.LCONE_ID theInputConeMosaic.MCONE_ID theInputConeMosaic.SCONE_ID];
    metaData.coneColors = [theInputConeMosaic.lConeColor; theInputConeMosaic.mConeColor; theInputConeMosaic.sConeColor];
    
    % Source lattice (i.e. cone mosaic lattice) struct
    sourceLatticeStruct = struct(...
        'name', 'cone RFs', ...
        'DegsToMMsConversionFunction', customDegsToMMsConversionFunction, ...
        'MMsToDegsConversionFunction', customMMsToDegsConversionFunction, ...
        'RFpositionsMicrons', theInputConeMosaic.coneRFpositionsMicrons, ...
        'metaData', metaData ...
        ); 
    

    % Import mRGC RF positions (destination) 
    mRGCRFposMicrons = retinalattice.import.finalMRGCPositions(...
                 sourceLatticeSizeDegs, ...
                 mean(sourceLatticeStruct.RFpositionsMicrons,1), ... 
                 max(max(sourceLatticeStruct.RFpositionsMicrons,[], 1)-min(sourceLatticeStruct.RFpositionsMicrons,[], 1)), ...
                 theInputConeMosaic.whichEye, ...
                 customDegsToMMsConversionFunction);


    % Destination lattice struct
    destinationLatticeStruct = struct(...
        'name', 'mRGC RFs', ...
        'DegsToMMsConversionFunction', customDegsToMMsConversionFunction, ...
        'MMsToDegsConversionFunction', customMMsToDegsConversionFunction, ...
        'RFpositionsMicrons', mRGCRFposMicrons ...
        );

    % Instantiate a coneToMidgetRGCConnector that will only connect
    % L and M-cones to RGCs
    lmConeIndices = [...
        theInputConeMosaic.lConeIndices(:); ...
        theInputConeMosaic.mConeIndices(:)];

    theMidgetRGCconnectorOBJ = coneToMidgetRGCConnector(...
        sourceLatticeStruct, destinationLatticeStruct, ...
        'coneIndicesToBeConnected', lmConeIndices, ...
        'visualizeConnectivityAtIntermediateStages', true);
  

end
