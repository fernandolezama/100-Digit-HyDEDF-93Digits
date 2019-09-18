% deParameters.I_NP= 20; % population in DE
% deParameters.F_weight= 0.3; %Mutation factor
% deParameters.F_CR= 0.5; %Recombination constant
% deParameters.I_itermax= 500; % number of max iterations/gen
% deParameters.I_strategy   = 1; %DE strategy
% deParameters.adaptActivated=1; %Change this to one for adaptive DE

deParameters.I_bnd_constr = 3; %Using bound constraints /is possible to change direct in DE
% 1 repair to the lower or upper violated bound (why)
% 2 rand value in the allowed range
% 3 bounce back

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Settings used in WCCI paper
deParameters.adaptActivated=0; %Change this to one for adaptive DE
deParameters.I_strategy=13;
deParameters.I_strategyVersion=3;
deParameters.I_itermax= 1e3;
            
deParameters.I_NP=50;
deParameters.F_weight=0.5;
deParameters.F_CR=0.9;

% deParameters.FirstEval=10; %First evaluation
% %deParameters.SecondEval=10; %Remaining evaluations
% 
% %deParameters.I_itermax= 100+2000; % number of max iterations/gen
% deParameters.I_itermax=floor(50000/(deParameters.I_NP*deParameters.FirstEval))-1;
% 
% 
% deParameters.Criterion=3; 
%1: Worst criterion
%2: Best Value
%3: mean value

%deParameters.fnc='fitnessFun_DER_WCCI';
%otherParameters = setOtherParameters(caseStudyData(1),deParameters.I_NP); %By FLC      
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%