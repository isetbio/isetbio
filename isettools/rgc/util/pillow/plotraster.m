function [MtspPlot, hh] = plotraster(varargin);
% function [MtspPlot, hh] = plotraster(Mtsp1, ... Mtspn, trange, rpts, lw);
% 
% Produces a raster plot of the spike time data contained in Mtsp.
% Inputs:  
%         Mtsp = matrix of spike times.  Each column should contain
%            the spike times for a single repeat of the stimulus (padded
%            with zeros if necessary).
% 
%            (Alternative:  Mtsp can be an array of 1's and 0's.  In
%            this case, each row of Mtsp is interpreted as one repeat of
%            the stimulus and time is taken in units of 1 per column of
%            Mtsp).
%
%         trange = [t0 t1], plot all spikes after time t0, before time
%            t1.  If trange a scalar, assumes t0 = 0.
%
%         rpts = [i1 in], which columns of Mtsp to use  (OPTIONAL)
%            If rpts a scalar, assumes i1 = 1.
%
%         lw = linewidth for plotted lines  (matlab default is 1).
%
% Outputs: 
%          MtspPlot = subset of spike times plotted for current figure
%          hh = handles to lines


lw0 = 1;  % default line width

%  Parse input to allow for multiple Mtsp arguments
maxn = 0;
for j = 1:nargin
    if max(size(varargin{j})) > 2
	k =j;
	if maxn < size(varargin{j},1)	    
	    maxnsp = size(varargin{j},1);
	end
	
    end
end

if k > 1
    Mtsp = zeros(maxnsp,size(varargin{j},1));
    j0 = 0;
    for j = 1:k
	[m,n] = size(varargin{j});
	Mtsp(1:m,j0+1:j0+n) = varargin{j};
	j0 = j0+n;
    end
else
    Mtsp = varargin{1};
end
varargin(1:k) = [];

if iscell(Mtsp)
    Mtsp = convrtTspCelltoMatrix(Mtsp);
end

% Set other pars
switch (nargin-k)
 case 0,
  trnge = [0 max(max(Mtsp))];
  irnge = [1 size(Mtsp,2)];
  lw = lw0;
 case 1,
  trnge = varargin{1};
  irnge = [1 size(Mtsp,2)];
  lw = lw0;
 case 2
  trnge = varargin{1};
  irnge = varargin{2};
  lw = lw0;
 case 3,
  trnge = varargin{1};
  irnge = varargin{2};
  lw = varargin{3};
 otherwise,
   fprintf(1, 'STOPPING--Check arguments to plotraster.m\n');
   keyboard;
end

if length(trnge) < 2
    trnge = [0 trnge];
end
if length(irnge) < 2
    irnge = [1 irnge];
end

% Check to make sure Mtsp is appropriate for rpts
if irnge(2) > size(Mtsp,2)
    error('Too few columns in Mtsp for requested range of repeats');
end

MtspPlot = zeros(size(Mtsp,1), irnge(2)-irnge(1)+1);
maxnsp = 0;

hldstate = 0;
if ishold
    hldstate = 1;
end

hh = [];
%  Create plot
for j = irnge(1):irnge(2)
    sptimes = Mtsp(find((Mtsp(:,j)>trnge(1)) & (Mtsp(:,j) <= trnge(2))),j);
    nsp = length(sptimes);
    ypl = repmat([j-irnge(1); j-irnge(1)+.9], 1,  nsp);   
    % hh(end+1:end+length(ypl)) = plot([sptimes'; sptimes'], ypl, 'k', 'linewidth', lw); hold on;
        
    MtspPlot(1:nsp,j-irnge(1)+1) = sptimes;
    if nsp > maxnsp
	maxnsp = nsp;
    end       
end
if ~hldstate
    hold off;
end

% % set(gca, 'xlim', [trnge(1) trnge(2)]); 
% axis([trnge(1) trnge(2) 0 irnge(2)-irnge(1)+1]);  axis ij;
% MtspPlot(maxnsp+1:end,:) = [];  % remove extra rows from MtspPlot
% set(gca, 'tickdir', 'out');

%  -------------- Cut text from above ---------------
% % test non-monotonicity of spike times in first column of Mtsp, to see
% % if we have to convert array to spike times
% m1 = Mtsp(:,1); 
% if (any(diff(m1(1:max(find(m1)))) <= 0))
%     fprintf(1, 'plotraster:  Converting spike array to matrix of spike times\n');
%     nsps = max(sum(Mtsp'));
%     m = size(Mtsp,1);
%     M2 = zeros(nsps, size(Mtsp,1));
%     for j = 1:m
% 	sptimes = find(Mtsp(j,:));
% 	M2(1:length(sptimes),j) = sptimes';
%     end    
%     Mtsp = M2;
% end
% 	
