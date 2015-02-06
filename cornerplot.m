function [fig,ax]=cornerplot(data, varargin)
%CORNERPLOT Corner plot showing projections of a multidimensional data set.
%
% CORNERPLOT(DATA) plots every 2D projection of a multidimensional data
% set. DATA is nSamples by nDimensions matrix.
%
% CORNERPLOT(DATA,NAMES) prints the names of each dimension. NAMES is a
% cell array of strings of length nDims, or an empty cell array.
%
% CORNERPLOT(DATA,NAMES,TRUTHS) indicates reference values on the plots.
% TRUTHS is a vector of length nDims. This might be useful, for instance,
% when looking at samples from fitting to synthetic data, where true
% parameter values are known.
%
% CORNERPLOT(DATA,NAMES,TRUTHS,BOUNDS) indicates lower and upper bounds for
% each dimension. BOUNDS is a 2 by nDims matrix where the first row is the
% lower bound for each dimension, and the second row is the upper bound.
%
% FIG is the handle for the figure, and AX is a
% nDimensions-by-nDimensions array of all subplot handles.
%
% Inspired by the more excellent and full-featured triangle.py
% (github.com/dfm/triangle.py) by Dan Foreman-Mackey (dan.iel.fm).
%
% Requires kde2d
% (mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation/content/kde2d.m)
% by Zdravko Botev (web.maths.unsw.edu.au/~zdravkobotev/).
%
% William Adler, January 2015
% Ver 1.0
% will@wtadler.com

if length(size(data))~=2
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
bounds_supplied = true;

if nargin > 1
    names = varargin{1};
    if ~isempty(names) && ~(iscell(names) && length(names)==nDims)
        error('NAMES must be a cell array with length equal to the number of dimensions in your data.')
    end
    if nargin > 2
        truths = varargin{2};
        if ~isempty(truths) && ~(isfloat(truths) && numel(truths)==nDims)
            error('TRUTHS must be a vector with length equal to the number of dimensions in your data.')
        end
        if nargin > 3
            bounds = varargin{3};

            if ~isempty(bounds) && ~(isfloat(bounds) && all(size(bounds)==[2 nDims]))
                error('BOUNDS must be a 2-by-nDims matrix.')
            end
        end
    end
end

if isempty(bounds)
    bounds = nan(nDims,2);
    bounds_supplied = false;
end

% plotting parameters
fig=figure;
ax=nan(nDims);
hist_bins = 100;
lines = 10;
res = 2^6; % defines grid for which kde2d will compute density. must be a power of 2.
linewidth = 1;
axes_defaults = struct('tickdirmode','manual',...
    'tickdir','out',...
    'ticklength',[.035 .035],...
    'box','off');


% plot histograms
for i = 1:nDims
    if ~bounds_supplied
        bounds(:,i) = [min(data(:,i)) max(data(:,i))];
    end
    
    ax(i,i)=subplot_tight(1+nDims,1+nDims, sub2ind([1+nDims 1+nDims],1+i,i));
    [n,x]=hist(data(:,i), hist_bins);
    plot(x,n/sum(n),'k-')
    set(gca,'yticklabel',[],'xlim',bounds(:,i),'ylim',[0 max(n/sum(n))],axes_defaults);
    
    if i~=nDims
        set(gca,'xticklabel',[])
    end
    
    if ~isempty(truths)
        hold on
        plot([truths(i) truths(i)], [0 1], 'k-', 'linewidth',linewidth)
    end
    
    if ~isempty(names)
        ylabel(sprintf('p(%s)',names{i}))
        if i==nDims
            xlabel(names{i})
        end
    end
    
end

% plot projections
if nDims > 1
    for p1 = 1:nDims-1
        for p2 = p1+1:nDims
            [~, density, X, Y] = kde2d([data(:,p1) data(:,p2)],res,[bounds(1,p1) bounds(1,p2)],[bounds(2,p1) bounds(2,p2)]);
            ax(p2,p1)=subplot_tight(1+nDims,1+nDims, sub2ind([1+nDims 1+nDims],1+p1,p2));
            contour(X,Y,density, lines)
            

            set(gca,'yticklabel',[],'xticklabel',[],'xlim',bounds(:,p1),'ylim',bounds(:,p2), axes_defaults);

            if ~isempty(truths)
                yl = get(gca,'ylim');
                xl = get(gca,'xlim');
                hold on
                plot(xl, [truths(p2) truths(p2)],'k-', 'linewidth',linewidth)
                plot([truths(p1) truths(p1)], yl,'k-', 'linewidth',linewidth)
            end
            if p1==1
                if ~isempty(names)
                    ylabel(names{p2})
                end
                set(gca,'yticklabelmode','auto')
            end
            if p2==nDims
                if ~isempty(names)
                    xlabel(names{p1})
                end
                set(gca,'xticklabelmode','auto')
            end
        end
        
        % link axes
        row = ax(1+p1,:);
        row = row(~isnan(row));
        row = row(1:end-1);
        
        col = ax(:,p1);
        col = col(~isnan(col));
        
        linkaxes(row, 'y');
        linkaxes(col, 'x');

    end
end
end

