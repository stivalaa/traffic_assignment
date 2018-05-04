function latency = BPR( x, xdata )
%BPR - Bureau of Public Roads latency function for use with lsqcurvefit
%   latency =  (1 + alpha * volratio^beta)
%
%   x      contains the parameters: x(1) = alpha, x(2) = beta
%   xdata is the vector of volume ratio (volratio values) Q/Qmax
%
% $Id: BPR.m 406 2011-06-23 05:12:53Z astivala $
%

latency = (1 + x(1) * xdata.^x(2));

end

