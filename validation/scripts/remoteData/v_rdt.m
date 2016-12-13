function varargout = v_rdt(varargin)
%
%  Validate ISETBIO-based colorimetric computations by comparing to PTB-based colorimetric computations.
%
    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

function ValidationFunction(runTimeParams)

%% Test remote data toolbox access

rd = RdtClient('isetbio');          % Open up the channel
rd.crp('/validation/full/color');   % Change to the relevant directory

a = rd.listArtifacts('print',true); % Maybe not needed.

assert(~isempty(a));                % There are a couple of things in there

foo = rd.readArtifact(a(1));

end

    

