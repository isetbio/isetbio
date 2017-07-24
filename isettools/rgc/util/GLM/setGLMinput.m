function glminput = setGLMinput(mosaic)
% DEPRECATED
%
% Reformats a mosaic's linear response into the simGLM or simGLMcpl format
%
% This is used for running Pillow's code.  At the moment we use it for
% other models, but that should be eliminated over time.
% 
% % \pillow_code_GLM_v1_Feb2010, available at 
% http://pillowlab.princeton.edu/code_GLM.html.
% 
% 
% 3/2016 JRG (c) isetbio

%% Set glminput to linear response

% nCells = mosaic.get('mosaic samples');
% nCellsTotal = nCells(1)*nCells(2);
% nSamples = size(lResponse,3);

lResponse = mosaic.get('response linear');
glminput = RGB2XWFormat(lResponse);

% glminput = zeros(nCellsTotal,nSamples));
% nanflag = 0; %#ok<*NASGU>

% cellCtr = 0;
% for i = 1:nCells(1)
%     for j = 1:nCells(2)        
%         cellCtr = cellCtr+1;
%         %  glminput(cellCtr,:) = (squeeze(sum(sum(response{i,j}(:,:,:,1),1),2)))';
%         
%         if sum(isnan(mosaic{i,j,1}))>0
%             nanflag = 1;
%             mosaic{i,j,1}(isnan(mosaic{i,j,1})) = 0;
%         end
%         
%         glminput(cellCtr,:) = mosaic{i,j,1};
%     end
% end

% if nanflag;
%     warning('NaN values in input replaced with zeros.');
% end

end
