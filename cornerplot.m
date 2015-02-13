function [fig,ax]=cornerplot(data, varargin)
%CORNERPLOT Corner plot showing projections of a multidimensional data set.
%
% CORNERPLOT(DATA) plots every 2D projection of a multidimensional data
% set. DATA is an nSamples-by-nDimensions matrix.
%
% CORNERPLOT(DATA,NAMES) prints the names of each dimension. NAMES is a
% cell array of strings of length nDimensions, or an empty cell array.
%
% CORNERPLOT(DATA,NAMES,TRUTHS) indicates reference values on the plots.
% TRUTHS is a vector of length nDimensions. This might be useful, for instance,
% when looking at samples from fitting to synthetic data, where true
% parameter values are known.
%
% CORNERPLOT(DATA,NAMES,TRUTHS,BOUNDS) indicates lower and upper bounds for
% each dimension. BOUNDS is a 2 by nDimensions matrix where the first row is the
% lower bound for each dimension, and the second row is the upper bound.
%
% FIG is the handle for the figure, and AX is a
% nDimensions-by-nDimensions array of all subplot handles.
%
% Inspired by triangle.py (github.com/dfm/triangle.py)
% by Dan Foreman-Mackey (dan.iel.fm).
%
% Requires kde2d
% (mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation/content/kde2d.m)
% by Zdravko Botev (web.maths.unsw.edu.au/~zdravkobotev/).
%
% William Adler, January 2015
% Ver 1.0
% will@wtadler.com

if ~exist('kde2d','file')
    error('You must install <a href="http://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation/content/kde2d.m">kde2d.m</a> by <a href="http://web.maths.unsw.edu.au/~zdravkobotev/">Zdravko Botev</a>.')
end

if length(size(data)) ~= 2
    error('x must be 2D.')
end

nDims = min(size(data));

% make sure columns are the dimensions of the data
if nDims ~= size(data,2)
    data = data';
end

% assign names and truths if given
names = {};
truths = [];
bounds = [];
bounds_supplied = true;

if nargin > 1
    names = varargin{1};
    if ~isempty(names) && ~(iscell(names) && length(names) == nDims)
        error('NAMES must be a cell array with length equal to the number of dimensions in your data.')
    end
    if nargin > 2
        truths = varargin{2};
        if ~isempty(truths) && ~(isfloat(truths) && numel(truths) == nDims)
            error('TRUTHS must be a vector with length equal to the number of dimensions in your data.')
        end
        if nargin > 3
            bounds = varargin{3};
            
            if ~isempty(bounds) && ~(isfloat(bounds) && all(size(bounds) == [2 nDims]))
                error('BOUNDS must be a 2-by-nDims matrix.')
            end
        end
    end
end

if isempty(bounds)
    bounds = nan(2,nDims);
    bounds_supplied = false;
end

% plotting parameters
fig = figure;
ax = nan(nDims);
hist_bins = 20;
lines = 10;
res = 2^6; % defines grid for which kde2d will compute density. must be a power of 2.
linewidth = 1;
axes_defaults = struct('tickdirmode','manual',...
    'tickdir','out',...
    'ticklength',[.035 .035],...
    'box','off',...
    'color','none');

% plot histograms
for i = 1:nDims
    if ~bounds_supplied
        bounds(:,i) = [min(data(:,i)) max(data(:,i))];
    end
    
    ax(i,i) = tight_subplot(1+nDims,1+nDims, i, 1+i);
    [n,x] = hist(data(:,i), hist_bins);
    plot(x,n/sum(n),'k-')
    set(gca,'yticklabel',[],'xlim',bounds(:,i),'ylim',[0 max(n/sum(n))],axes_defaults);
    
    if i ~= nDims
        set(gca,'xticklabel',[])
    end
    
    if ~isempty(truths)
        hold on
        plot([truths(i) truths(i)], [0 1], 'k-', 'linewidth',linewidth)
    end
    
    if ~isempty(names)
        if i == 1
            ylabel(names{i});
        end
        if i == nDims
            xlabel(names{i})
        end
    end
    
end

% plot projections
if nDims > 1
    for d1 = 1:nDims-1 % col
        for d2 = d1+1:nDims % row
            [~, density, X, Y] = kde2d([data(:,d1) data(:,d2)],res,[bounds(1,d1) bounds(1,d2)],[bounds(2,d1) bounds(2,d2)]);

            ax(d2,d1) = tight_subplot(1+nDims,1+nDims, d2, 1+d1);
            contour(X,Y,density, lines)
            
            set(gca,'yticklabel',[],'xticklabel',[],'xlim',bounds(:,d1),'ylim',bounds(:,d2), axes_defaults);
            
            if ~isempty(truths)
                yl = get(gca,'ylim');
                xl = get(gca,'xlim');
                hold on
                plot(xl, [truths(d2) truths(d2)],'k-', 'linewidth',linewidth)
                plot([truths(d1) truths(d1)], yl,'k-', 'linewidth',linewidth)
            end
            if d1 == 1
                if ~isempty(names)
                    ylabel(names{d2})
                end
                set(gca,'yticklabelmode','auto')
            end
            if d2 == nDims
                if ~isempty(names)
                    xlabel(names{d1})
                end
                set(gca,'xticklabelmode','auto')
            end
        end
        
        % link axes
        row = ax(1+d1,:);
        row = row(~isnan(row));
        row = row(1:end-1);
        
        col = ax(:,d1);
        col = col(~isnan(col));
        
        linkaxes(row, 'y');
        linkaxes(col, 'x');
        
    end
end
end

function h=tight_subplot(m,n,subplot_row,subplot_col)
% adapted from subplot_tight
% (mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot/)
% by Nikolay S. (http://vision.technion.ac.il/~kolian1/)

margins = [.015, .015];

height = (1-(m+1)*margins(1))/m;
width = (1-(n+1)*margins(2))/n;

bottom =(m-subplot_row)*(height+margins(1))+margins(1); % merged subplot bottom position
left =subplot_col*(width+margins(2))-width;              % merged subplot left position
pos_vec=[left bottom width height];

h=subplot('Position',pos_vec);
end