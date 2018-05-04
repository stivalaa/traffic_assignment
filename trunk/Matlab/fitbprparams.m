%
% File:    fitbprparams.m
% Author:  Alex Stivala
% Created: June 2011
%
%
% MATLAB script to fit BPR function to SpeedFlowVoyager data.
% This data is a table with some data points giving speed (km/h)
% for different volume rations (volume/capacity) from 0 to 1.1 (10% over
% capacity)
% We convert these to latency and fit alpha, beta
% parameters to the BPR function
% We converted the SpeedflowVoyager matrix from CSV in such a way that
% the first row is just the index (0 for first col, capindex 2 is col
% 2, etc.), so not really used. The first colum is the volratio
% and each subsequent colum is colum vector where each element is
% speed corresponing to volratio in that row.
%
% Writes the output file bprparams.matrix in ASCII MATLAB matrix format
% WARNING: overwrites if it exsits.
%
% $Id: fitbprparams.m 407 2011-06-23 05:31:38Z astivala $
%

load SpeedflowVoyager;
[nrows, ncols] = size(SpeedflowVoyager);
volratios = SpeedflowVoyager(2:nrows, 1);
parammatrix = zeros(3, ncols); % each col has 2 params+freeflow speed for that cols index
for capindex = 2:ncols
    x0 = [0.15; 4.0]; % 'standard' paranmeters
    options = optimset('Display','off');
    params = lsqcurvefit(@BPR, x0, volratios,                          ...
    SpeedflowVoyager(2,capindex)./SpeedflowVoyager(2:nrows,capindex),  ...
                         [], [], options);
    parammatrix(:,capindex) = [params; SpeedflowVoyager(2, capindex)];
    
    figure(capindex);
    plot(volratios, SpeedflowVoyager(2,capindex)  ./                    ... 
                    SpeedflowVoyager(2:nrows,capindex), 'b');
    hold;
    xp = 0:0.01:1.2;
    plot(xp, BPR(parammatrix(:,capindex), xp), 'r');
    
end;

% Actually seems better to take the transpose of parammatrix so that
% we have only 3 columns (the 2 parameters + max speed) and 1 row for each capindex
% (note first row is all 0, capindex starts at 2)

bprparams = parammatrix';
save('bprparams.matrix','bprparams','-ascii')
