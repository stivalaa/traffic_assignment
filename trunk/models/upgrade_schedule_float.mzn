%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:    upgrade_schedule_float.mzn
% Author:  ADS
% Created: July 2011
%
%
% $Id: upgrade_schedule_float.mzn 614 2011-08-25 00:15:21Z astivala $
%
% MiniZinc model to find optimum schedule of road upgrades given, for each
% year, the
% cost of each upgrade, and its benefit in terms of vehicle hours
% travelled (VHT), as well as the difference between the sum of each
% pair of upgrades benefits and the actual benefit of the pair
% which is estimated as 0 (i.e. pair of upgrades gives benefit of sum
% of each separate upgrade) except where we have deteremined that this
% pair is not disjoint so we have some difference between the sum of
% individual benefits and the benefit of the pair.
% this data is obtained from the output of "traffic assignment problem"
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

% difference in pairwise versus sum of indiviudal VHT reduction,
% in each period
array[Periods, Upgrades, Upgrades] of float: PairMinusSumVHT =
    array3d(Periods, Upgrades, Upgrades, PairMinusSumVHT_list);
array[1..n*n*numperiods] of float: PairMinusSumVHT_list;



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

% The objective function is quadratic because as well as the benefit
% for each individual upgrade we account for the possibly different
% benefit for each pair of upgrades in the case both in the pair are in
% the subset considered



var float: totalBenefitDollarsNPV = 
    sum(t in Periods)(
        sum (u,v in Upgrades where u < v)(
            (int2float(bool2int(u in upgradeSet[t])) *
             int2float(bool2int(v in upgradeSet[t])) *
            VHT2dollarFactor *
            PairMinusSumVHT[t,u,v] / pow(1.0+interestRate,int2float(t))
        ) ) 
        +
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
