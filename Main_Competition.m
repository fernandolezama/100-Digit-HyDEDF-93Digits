clear;close all;clc;

%% This package is a MATLAB source code of HyDE-DF with reinitialization strategy.
%% About HyDE-DF 1.2, please see following papers:
%%
%% * Lezama et. al: Hybrid-adaptive differential evolution with decay function (HyDE-DF) applied to the 100-digit challenge competition on single objective numerical optimization. In Proceedings of the Genetic and Evolutionary Computation Conference Companion (GECCO '19). 2019 DOI: https://doi.org/10.1145/3319619.3326747
%%
%% Notify that line this HyDE-DF 1.2 implement a different strategy for reinitialization compared with HyDE-DF to be published
%% Contact ing.flezama@gmail.com for details in the implementation

rand('seed', sum(100 * clock));
for alg_test=1:1

addpath('HyDE')
% Algorithm parameters
DEparameters

if alg_test==1 %This is the best configuration found so far
    deParameters.I_strategy=3;
    deParameters.I_strategyVersion=2;
    deParameters.I_itermax_DE=1e5; %Duration of vortex function
    
    deParameters.I_itermax= 1e7; %Maximum budget

    deParameters.I_NP=50;
    deParameters.F_weight=0.5;
    deParameters.F_CR=0.5;
    deParameters.I_bnd_constr=3;
end


FUNCanalysis=9;
for j=9:FUNCanalysis

 FN=j
 if FN==7 || FN==9
     deParameters.I_NP=1000;
 end
 if FN==8 
    deParameters.I_NP=200;
 end
 if FN==10
     deParameters.I_NP=50;
 end
 
%Function parameters
func = callFunction(FN); %get the function struct for [1...50] functions
otherParameters.objfun = func.name; % function to be optimized
otherParameters.objfunCode=func.code;
otherParameters.dim = func.dim; %dimension of the problem
otherParameters.lowerlimit =  func.lowerlimit; %lower limit of the problem
otherParameters.upperlimit = func.upperlimit; %upper limit of the problem
lowerlimit =  func.lowerlimit; %lower limit of the problem
upperlimit = func.upperlimit; %upper limit of the problem

tTotalTime=tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set lower/upper bounds of variables 
%[lowerB,upperB] = setVariablesBounds(caseStudyData,otherParameters);

lowerB=lowerlimit*ones(1,otherParameters.dim);
upperB=upperlimit *ones(1,otherParameters.dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call the MH for optimizationclear 
noRuns = 50; % Number of trials here
ResDB=struc([]);
%M_Enhance=struc([]);
   parfor iRuns=1:noRuns %Number of trails
        tOpt=tic;
        rand('state',sum(iRuns*100*clock))% ensure stochastic indpt trials
           [ResDB(iRuns).Fit_and_p, ...
              ResDB(iRuns).sol, ...
              ResDB(iRuns).fitVector,...
              ResDB(iRuns).table1,...
              ResDB(iRuns).enhance] =...
              HyDE(deParameters,otherParameters,lowerB,upperB);
          

       ResDB(iRuns).tOpt=toc(tOpt); % time of each trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save the results and stats
%        Save_results
    iRuns
    end
tTotalTime=toc(tTotalTime); %Total time


%General stats
for i=1:noRuns
    Values(i)=  ResDB(i).Fit_and_p;
end
alf=1:3:3*50;
Summary(alf(FN):alf(FN)+2,1)=[mean(Values);std(Values);min(Values) ];

for i=1:noRuns
    Table1(:,i,FN)=  ResDB(i).table1;
end

Count=sum(Table1(1:10,:,FN)~=0);

Score=sort(Count,'descend');
Score_f=sum(Score(1:25))/25;

for i=1:11
    Table2(FN,i)=sum(Count==i-1);
end
Table2(FN,12)=Score_f;


if alg_test==1
    filename=['Results_HyDEDF_S00/funct_'  num2str(FN)];
end

save(filename,'-v7.3')

%% End of MH Optimization
end

if alg_test==1
  save('Results_HyDEDF_S00/TableX','Summary')
  save('Results_HyDEDF_S00/Table1','Table1')
  save('Results_HyDEDF_S00/Table2','Table2')
end


end

%addpath('Functions') 
% CEC19 Test Function Suite 
%   Noor Awad (email: noorawad1989@gmail.com) 
%   Dec. 13th 2018
%   1. Run the following command in Matlab window:
%   mex cec19_func.cpp -DWINDOWS
%   2. Then you can use the test functions as the following example:
%   f = cec19_func(x,func_num); 
%   Here x is a D*pop_size matrix.

