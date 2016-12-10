%% Test remote data toolbox access

rd = RdtClient('isetbio');
rd.crp('/validation/full/color');
a = rd.listArtifacts('print',true);
assert(length(a) == 2)

%%
    

