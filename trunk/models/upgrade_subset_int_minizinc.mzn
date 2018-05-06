%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:    upgrade_subset_int.zinc
% Author:  ADS
% Created: May 2011
%
%
% $Id: upgrade_subset_int_minizinc.mzn 347 2011-06-06 06:47:13Z astivala $
%
% Zinc model to find optimum subset of road upgrades given the
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
% Cannot use MIP solver to solve (get segfault if mzn-g12mip is run on
% this model), also mzn-g12lazy doesn't work (wrong answers) and error
% with mzn-g12sat. So must use mzn-g12fd (minizinc default?)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% data defined in .dzn file
%

int: n;                                % number of individual upgrades
int: Budget;                           % total budget
set of int: Upgrades = 1..n; 
array[Upgrades] of int: cost ;         % cost of each upgrade
array[Upgrades] of int: benefitVHT;    % reduction in VHT for each upgrade
array[Upgrades] of string: upgradeName;% name describing each upgrade
array[Upgrades,Upgrades] of int: PairMinusSumVHT; % difference in pairwise
                                                  % versus sum of indiviudal


%
% variables
%

var set of Upgrades: upgradeSet;  % the optimal subset of all upgrades

%
% constraints
%

constraint sum (u in Upgrades) 
               (bool2int(u in upgradeSet) * cost[u]) <= Budget;


%
% objective function
%

% The objective function is quadratic because as well as the benefit
% for each individual upgrade we account for the possibly different 
% benefit for each pair of upgrades in the case both in the pair are in
% the subset considered

var int: totalBenefitVHT = sum (u,v in Upgrades where u < v)
      ( bool2int(u in upgradeSet) * bool2int(v in upgradeSet) *
          PairMinusSumVHT[u,v] ) + 
      sum(u in Upgrades)(bool2int(u in upgradeSet) * benefitVHT[u]);


solve maximize totalBenefitVHT;

%
% output
% 

output ["totalBenefitVHT = ", show(totalBenefitVHT), "\n",
        show([upgradeName[i] | i in Upgrades where i in fix(upgradeSet)]),
         "\n",
         "totalCost = ", show( sum(i in Upgrades where i in fix(upgradeSet))
                               (cost[i]) ), 
         "\n"       
         ];




