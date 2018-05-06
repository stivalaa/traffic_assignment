%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:    upgrade_subset_degree3_float.zinc
% Author:  ADS
% Created: May 2011
%
%
% $Id: upgrade_subset_degree3_float.mzn 637 2011-08-31 04:49:57Z astivala $
%
% Zinc model to find optimum subset of road upgrades given the
% cost of each upgrade, and its benefit in terms of vehicle hours
% travelled (VHT), as well as the difference between the sum of each
% pair of upgrades benefits and the actual benefit of the pair,
% AND the difference between sum of each triple of upgrades and
% actual benefit of the triple. So the objective function here
% is cubic.
% This data is obtained from the output of "traffic assignment problem"
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
% difference in pairwise versus sum of indiviudal
array[Upgrades,Upgrades] of float: PairMinusSumVHT;
% difference in 3-wise versus sum of individual and adjustments for pairs
% (due to PairMinusSumVHT). So maybe the name is a bit misleading since
% this value is not triple-[sum of indiviudal] but triple - 
% ([sum of individual] +  [sum of all PairMinusSumVHT for each included pair])
array[Upgrades,Upgrades,Upgrades] of float: TripleMinusSumVHT = 
    array3d(Upgrades,Upgrades,Upgrades,TripleMinusSumVHT_list);
array[1..n*n*n] of float: TripleMinusSumVHT_list;

%
% variables
%

var set of Upgrades: upgradeSet;  % the optimal subset of all upgrades

%
% constraints
%

constraint sum (u in Upgrades) 
               (int2float(bool2int(u in upgradeSet)) * cost[u]) <= Budget;


%
% objective function
%

% The objective function is cubic because as well as the benefit
% for each individual upgrade we account for the possibly different 
% benefit for each pair of upgrades in the case both in the pair are in
% the subset considered (quadratic term) and also for each triple of
% upgrades when all 3 in subset considered (cubic term)

var float: totalBenefitVHT = sum (u,v,w in Upgrades where u < v /\ v < w)
      (int2float(bool2int(u in upgradeSet)) *
       int2float(bool2int(v in upgradeSet)) *
       int2float(bool2int(w in upgradeSet)) *
          TripleMinusSumVHT[u,v,w] ) + 
      sum (u,v in Upgrades where u < v) 
      (int2float(bool2int(u in upgradeSet)) *
       int2float(bool2int(v in upgradeSet)) *
          PairMinusSumVHT[u,v] ) + 
      sum(u in Upgrades)(int2float(bool2int(u in upgradeSet)) * benefitVHT[u]);


solve maximize totalBenefitVHT;

% Use MIP solver to solve.
%solve::backend_mip(default) maximize totalBenefitVHT;


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


