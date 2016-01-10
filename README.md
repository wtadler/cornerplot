# cornerplot

cornerplot takes a nSamples-by-nDimensions array, and makes density plots of every combination of the dimensions, and histograms for each dimension.
This is especially useful when using MCMC; you can see how the parameters in your model interact, and whether there are any tradeoffs between them.

Inspired by [corner.py](https://github.com/dfm/corner.py) by [Dan Foreman-Mackey](http://dan.iel.fm/).


# Usage
<img src="http://wtadler.com/picdrop/cornerplot.png" width=30% align="left" />
You can obtain a plot like that pictured to the left with `cornerplot(randn(500, 3))`.
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />


<img src="http://wtadler.com/picdrop/cornerplot_labels.png" width=30% align="left" />

You can also label parameters and mark true parameter values, as in `cornerplot(randn(500, 3), {'a', 'b', 'c'}, randn(1, 3))`
<br />
<br />
<br />
<br />
<br />
<br />

# Installation
Just put cornerplot.m in your MATLAB path. You also need to install the free kernel density estimator [kde2d](http://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation) by [Zdravko Botev](http://web.maths.unsw.edu.au/~zdravkobotev/).
