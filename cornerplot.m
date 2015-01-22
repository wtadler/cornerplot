function [fig,ax]=cornerplot(data, varargin)
%CORNERPLOT Corner plot showing projections of a multidimensional data set.
% 
% CORNERPLOT(DATA) plots every 2D projection of a multidimensional data
% set. DATA is an array of size [nSamples, nDimensions].
% 
% CORNERPLOT(DATA,NAMES) prints the names of each dimension. NAMES is a
% cell array of strings of length nDims, or an empty cell array.
% 
% CORNERPLOT(DATA,NAMES,TRUTHS) indicates reference values on the plots.
% TRUTHS is a vector of length nDims. This might be useful, for instance,
% when looking at samples from fitting to synthetic data, where true
% parameter values are known.
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
nSamples=max(size(data));

% make sure columns are the dimensions of the data
if nDims ~= size(data,2)
    data = data';
end

% assign names and truths if given
names = {};
truths = [];
if nargin > 1
    names = varargin{1};
    if ~isempty(names) && ~(iscell(names) && length(names)==nDims)
        error('NAMES must be a cell array with length equal to the number of dimensions in your data.')
    end
    if nargin > 2
        truths = varargin{2};
        if ~(isfloat(truths) && numel(truths)==nDims)
            error('TRUTHS must be a vector with length equal to the number of dimensions in your data.')
        end
    end
end

% plotting parameters
fig=figure;
ax=nan(nDims);
hist_bins = 50;
lines = 10;
res = 2^6; % defines grid for which kde2d will compute density. must be a power of 2.
linewidth = 1;
axes_defaults = struct('tickdirmode','manual',...
    'tickdir','out',...
    'ticklength',[.035 .035],...
    'box','off');

lims = nan(nDims,2);

% plot histograms
for i = 1:nDims
    lims(i,:) = [min(data(:,i)) max(data(:,i))];
    
    ax(i,i)=subplot_tight(1+nDims,1+nDims, sub2ind([1+nDims 1+nDims],1+i,i));
    [n,x]=hist(data(:,i), hist_bins);
    plot(x,n/sum(n),'k-')
    set(gca,'yticklabel',[],'xlim',lims(i,:),axes_defaults);

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
            [~, density, X, Y] = kde2d([data(:,p1) data(:,p2)],res,[lims(p1,1) lims(p2,1)],[lims(p1,2) lims(p2,2)]);
            ax(p2,p1)=subplot_tight(1+nDims,1+nDims, sub2ind([1+nDims 1+nDims],1+p1,p2));
            contour(X,Y,density, lines)
            
            set(gca,'yticklabel',[],'xticklabel',[],'xlim',lims(p1,:),'ylim',lims(p2,:), axes_defaults);
            
            if ~isempty(truths)
                yl = get(gca,'ylim');
                xl = get(gca,'xlim');
                hold on
                plot(xl, [truths(p2) truths(p2)],'k-', 'linewidth',linewidth)
                plot([truths(p1) truths(p1)], yl,'k-', 'linewidth',linewidth)
            end
            if ~isempty(names)
                if p1==1
                    ylabel(names{p2})
                    set(gca,'yticklabelmode','auto')
                end
                if p2==nDims
                    xlabel(names{p1})
                    set(gca,'xticklabelmode','auto')
                end
            end
        end
    end
end