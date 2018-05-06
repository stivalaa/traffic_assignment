%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:    upgrade_schedule_independent.mzn
% Author:  ADS
% Created: July 2011
%
%
% $Id: upgrade_schedule_independent.mzn 716 2011-09-19 07:15:54Z astivala $
%
% MiniZinc model to find optimum schedule of road upgrades given, for each
% year, the
% cost of each upgrade, and its benefit in terms of vehicle hours
% travelled (VHT). In this model we assume all upgrades are indpendent.
% This data is obtained from the output of "traffic assignment problem"
% (TAP) from e.g. the tap_frankwolfe_mpi program in ../tap/.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% data defined in .dzn file
%

float:interestRate;                    % for present value calculations
float:VHT2dollarFactor;                % multiplier to convert VHT to dollars
int: n;                                % number of individual upgrades
int: numperiods;                       % number of time periods
array[1..numperiods] of float: Budget;   % total budget in each time period
                                         % (as present value)
set of int: Periods = 1..numperiods;

set of int: Upgrades = 1..n;
array[Upgrades] of float: cost ;         % cost of each upgrade (present value)
array[Periods,Upgrades] of float: benefitVHT;  % reduction in VHT for
                                                      % each upgrade in each
                                                      % time period
array[Upgrades] of string: upgradeName;% name describing each upgrade



%
% variables
%

% the optimal subset of available upgrades in each period
array[Periods] of var set of Upgrades: upgradeSet;


%
% constraints
%

% In each time period, we must fit within that period's total budget
constraint forall(t in Periods)
               (sum (u in Upgrades)
                (int2float( bool2int(u in upgradeSet[t]) )*cost[u]) <= Budget[t]
               ) ;


% Each upgrade can only be built once
constraint forall(u in Upgrades)(sum(t in Periods)
                                 (bool2int(u in upgradeSet[t])) <= 1);



%
% objective function
%



var float: totalBenefitDollarsNPV = 
    sum(t in Periods)(
        sum(u in Upgrades)(
            int2float(bool2int(u in upgradeSet[t])) *
            VHT2dollarFactor *
            benefitVHT[t,u] / pow(1.0+interestRate, int2float(t))
        )
    );

solve maximize totalBenefitDollarsNPV;



%
% output
%

output
    [ "totalBenefitDollarsNPV = ", show(totalBenefitDollarsNPV), "\n" ]
    ++
    [   "t = " ++ show(t) ++ ": " ++
        show([upgradeName[i] | i in Upgrades where i in fix(upgradeSet[t])])
        ++ "\t" ++
        "totalCost = " ++
        show(sum(i in Upgrades where i in fix(upgradeSet[t]))(cost[i])) ++
        "\n"
      | t in Periods ]
    ++
    [ "\n" ];
