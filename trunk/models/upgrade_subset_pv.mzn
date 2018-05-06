%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:    upgrade_subset_float.zinc
% Author:  ADS
% Created: May 2011
%
%
% $Id: upgrade_subset_pv.mzn 707 2011-09-15 06:48:33Z astivala $
%
% Zinc model to find optimum subset of road upgrades given the
% cost of each upgrade, and its benefit in terms of
% the dollar present value of vehicle hours
% travelled (VHT), as well as the difference between the sum of each
% pair of upgrades benefits and the actual benefit of the pair
% which is estimated as 0 (i.e. pair of upgrades gives benefit of sum
% of each separate upgrade) except where we have deteremined that this
% pair is not disjoint so we have some difference between the sum of
% individual benefits and the benefit of the pair.
% this data is obtained from the output of "traffic assignment problem"
% (TAP) from e.g. the tap_frankwolfe_mpi program in ../tap/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% data defined in .dzn file
%

int: n;                                % number of individual upgrades
float: Budget;                           % total budget
set of int: Upgrades = 1..n; 
array[Upgrades] of float: cost ;         % cost of each upgrade
array[Upgrades] of float: benefitVHT;    % reduction in VHT for each upgrade
array[Upgrades] of string: upgradeName;% name describing each upgrade
array[Upgrades,Upgrades] of float: PairMinusSumVHT; % difference in pairwise
                                                  % versus sum of indiviudal
float:interestRate;                    % for present value calculations
float:VHT2dollarFactor;                % multiplier to convert VHT to dollars
int:t;                                 % time period
set of Upgrades: potentialUpgradeSet; % subset of all upgrades to consider building
set of Upgrades: previouslyBuiltSet; % upgrades already built - account for them
                                         % computing VHT pairwise adjustments
                                         % but not in budget

%
% variables
%

var set of Upgrades: upgradeSet;  % the optimal subset of all upgrades

%
% constraints
%

% 
constraint sum (u in Upgrades) 
               (int2float(bool2int(u in upgradeSet)) * cost[u]) <= Budget;

constraint upgradeSet subset potentialUpgradeSet;

constraint assert( potentialUpgradeSet intersect previouslyBuiltSet == {},
           "nonempty intersection of potential and already built sets"  );

%
% objective function
%

% The objective function is quadratic because as well as the benefit
% for each individual upgrade we account for the possibly different 
% benefit for each pair of upgrades in the case both in the pair are in
% the subset considered


var float: totalBenefitDollarsPV = sum (u,v in Upgrades where u < v)
      (int2float(bool2int(u in upgradeSet /\ v in previouslyBuiltSet  \/
                          v in upgradeSet /\ u in previouslyBuiltSet  \/
                          u in upgradeSet /\ v in upgradeSet)) *
          PairMinusSumVHT[u,v] * 
          VHT2dollarFactor / pow(1.0+interestRate, int2float(t)) ) + 
      sum(u in Upgrades)(int2float(bool2int(u in upgradeSet)) * benefitVHT[u]
          * VHT2dollarFactor / pow(1.0+interestRate, int2float(t)));


solve maximize totalBenefitDollarsPV;



%
% output
% 

output ["totalBenefitDollarsPV = ", show(totalBenefitDollarsPV), "\n",
        show([upgradeName[i] | i in Upgrades where i in fix(upgradeSet)]),
         "\n",
         "totalCost = ", show( sum(i in Upgrades where i in fix(upgradeSet))
                               (cost[i]) ), 
         "\n"       
         ];


