%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB source code of HyDE-DF with reinitialization strategy.
%% About HyDE-DF 1.2, please see following papers:
%%
%% * Lezama et. al: Hybrid-adaptive differential evolution with decay function (HyDE-DF) applied to the 100-digit challenge competition on single objective numerical optimization. In Proceedings of the Genetic and Evolutionary Computation Conference Companion (GECCO '19). 2019 DOI: https://doi.org/10.1145/3319619.3326747
%%
%% Notify that line of code 220 is what makes HyDE-DF 1.2 different from HyDE-DF 1.1 to be published
%% Contact ing.flezama@gmail.com for details in the implementation


%Modification of DE
% Function:         [FVr_bestmem,S_bestval,I_nfeval] = deopt(fname,S_struct)
% Author:           Rainer Storn, Ken Price, Arnold Neumaier, Jim Van Zandt
% Modified by FLC \GECAD 2019

function [Fit_and_p,FVr_bestmemit, fitMaxVector,outcome,Enhance] = ...
    HyDE(deParameters,otherParameters,low_habitat_limit,up_habitat_limit,initialSolution)

outcome = zeros(11,1);
flat_Success=0;
digit_1_reach =  10^(0); digit_2_reach =  10^(-1); digit_3_reach =  10^(-2); digit_4_reach =  10^(-3); digit_5_reach =  10^(-4);
digit_6_reach =  10^(-5); digit_7_reach =  10^(-6); digit_8_reach =  10^(-7); digit_9_reach =  10^(-8); digit_10_reach =  10^(-9);

digit_1 = 0; digit_2 = 0; digit_3 = 0; digit_4 = 0; digit_5 = 0;
digit_6 = 0; digit_7 = 0; digit_8 = 0; digit_9 = 0; digit_10 = 0;


%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP;
F_weight     = deParameters.F_weight;
F_CR         = deParameters.F_CR;
I_D          = numel(up_habitat_limit); %Number of variables or dimension
deParameters.nVariables=I_D;
FVr_minbound = low_habitat_limit;
FVr_maxbound = up_habitat_limit;
I_itermax    = deParameters.I_itermax;

I_itermax_DE = deParameters.I_itermax_DE;
%reset_V=I_itermax_DE;

%Repair boundary method employed
BRM=deParameters.I_bnd_constr; %1: bring the value to bound violated
                               %2: repair in the allowed range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_strategy   = deParameters.I_strategy; %important variable
%fnc= otherParameters.objfun;
fnc= otherParameters.objfunCode;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----Check input variables---------------------------------------------
if (I_NP < 5)
   I_NP=5;
   fprintf(1,' I_NP increased to minimal value 5\n');
end
if ((F_CR < 0) || (F_CR > 1))
   F_CR=0.5;
   fprintf(1,'F_CR should be from interval [0,1]; set to default value 0.5\n');
end
if (I_itermax <= 0)
   I_itermax = 200;
   fprintf(1,'I_itermax should be > 0; set to default value 200\n');
end

%-----Initialize population and some arrays-------------------------------
%FM_pop = zeros(I_NP,I_D); %initialize FM_pop to gain speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-allocation of loop variables
fitMaxVector = nan(1,I_itermax);
% limit iterations by threshold
gen = 1; %iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----FM_pop is a matrix of size I_NPx(I_D+1). It will be initialized------
%----with random values between the min and max values of the-------------
%----parameters-----------------------------------------------------------
% FLC modification - vectorization
minPositionsMatrix=repmat(FVr_minbound,I_NP,1);
maxPositionsMatrix=repmat(FVr_maxbound,I_NP,1);
deParameters.minPositionsMatrix=minPositionsMatrix;
deParameters.maxPositionsMatrix=maxPositionsMatrix;

% generate initial population.
FM_pop=genpop(I_NP,I_D,minPositionsMatrix,maxPositionsMatrix);


if nargin>5
    noInitialSolutions = size(initialSolution,1);
    FM_pop(1:noInitialSolutions,:)=initialSolution;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Evaluate the best member after initialization----------------------
% Modified by FLC
f = cec19_func(FM_pop',fnc);
OFE=I_NP;
S_val=f';

[~,I_best_index] = min(S_val); % This mean that the best individual correspond to the best worst performance
FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
fitMaxVector(:,1)=S_val(I_best_index); %We save the mean value and mean penalty value
%------DE-Minimization---------------------------------------------
%------FM_popold is the population which has to compete. It is--------
%------static through one iteration. FM_pop is the newly--------------
%------emerging population.----------------------------------------
FVr_rot  = (0:1:I_NP-1);               % rotating index array (size I_NP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HYDE
if deParameters.I_strategy==3;
        F_weight_old=repmat(F_weight,I_NP,3);
        F_weight= F_weight_old;
        F_CR_old=repmat(F_CR,I_NP,1);
        F_CR=F_CR_old;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_strategyVersion=deParameters.I_strategyVersion;

%other.a=(I_itermax_DE-gen)/I_itermax_DE;
other.a_vector=zeros(1,I_itermax);
other.a_vector(1:I_itermax_DE)=(I_itermax_DE-(1:I_itermax_DE))./I_itermax_DE;
%other.a_vector(1:I_itermax)=(I_itermax-(1:I_itermax))./I_itermax;


cont=1;
other.lowerlimit=otherParameters.lowerlimit; %lower limit of the problem
other.upperlimit = otherParameters.upperlimit; %upper limit of the problem

while gen<=I_itermax  %%&&  fitIterationGap >= threshold
    
    %a = itr / MaxItr; % a value for gammaincinv function
    other.a=other.a_vector(cont); %Select value for decreasing functions
    
    
     if deParameters.I_strategy==3;
        value_R=rand(I_NP,3);
        ind1=value_R<0.1;
        ind2=rand(I_NP,1)<0.1;
        F_weight(ind1)=0.1+rand(sum(sum(ind1)),1)*0.9;
        F_weight(~ind1)=F_weight_old(~ind1);
        F_CR(ind2)=rand(sum(ind2),1);
        F_CR(~ind2)=F_CR_old(~ind2);
     end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [FM_ui,FM_base,~]=generate_trial(I_strategy,F_weight, F_CR, FM_pop, FVr_bestmemit,I_NP, I_D, FVr_rot,I_strategyVersion,other);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Boundary Control
    FM_ui=update(FM_ui,minPositionsMatrix,maxPositionsMatrix,BRM,FM_base);

    %Evaluation of new Pop
    %S_val_temp=feval(fnc,FM_ui);
    f = cec19_func(FM_ui',fnc);
    OFE=OFE+I_NP;%+I_NP;
    S_val_temp=f';
    
    
    %% Elitist Selection
    ind=S_val_temp<=S_val;
    S_val(ind)=S_val_temp(ind);
    FM_pop(ind,:)=FM_ui(ind,:);
  
    
    %% update best results
    [S_bestval,I_best_index] = min(S_val);
    FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
    % store fitness evolution and obj fun evolution as well
    fitMaxVector(1,gen) = S_bestval;
    %S_bestval
    
     if deParameters.I_strategy==3; %jDE
        F_weight_old(ind,:)=F_weight(ind,:);
        F_CR_old(ind)=F_CR(ind);
     end
    
     bsf_error_val = S_bestval - 1;
     
        %% record digits and FEs
        if bsf_error_val <= digit_1_reach  && digit_1 == 0
            fprintf('1st %d Gen: %d\n', OFE,gen);digit_1 = 1;outcome(1)=OFE;end
        if bsf_error_val <= digit_2_reach  && digit_2 == 0
            fprintf('2nd %d Gen: %d\n', OFE,gen);digit_2 = 1;outcome(2)=OFE;end
         if bsf_error_val <= digit_3_reach  && digit_3 == 0
            fprintf('3th %d Gen: %d\n', OFE,gen);digit_3 = 1;outcome(3)=OFE;end
         if bsf_error_val <= digit_4_reach  && digit_4 == 0
            fprintf('4th %d Gen: %d\n', OFE,gen);digit_4 = 1;outcome(4)=OFE;end
         if bsf_error_val <= digit_5_reach  && digit_5 == 0
            fprintf('5th %d Gen: %d\n', OFE,gen);digit_5 = 1;outcome(5)=OFE;end
         if bsf_error_val <= digit_6_reach  && digit_6 == 0
            fprintf('6th %d Gen: %d\n', OFE,gen);digit_6 = 1;outcome(6)=OFE;end
        if bsf_error_val <= digit_7_reach  && digit_7 == 0
            fprintf('7th %d Gen: %d\n', OFE,gen);digit_7 = 1;outcome(7)=OFE;end
         if bsf_error_val <= digit_8_reach  && digit_8 == 0
            fprintf('8th %d Gen: %d\n', OFE,gen);digit_8 = 1;outcome(8)=OFE;end
         if bsf_error_val <= digit_9_reach  && digit_9 == 0
            fprintf('9th %d Gen: %d\n', OFE,gen);digit_9 = 1;outcome(9)=OFE;end
         if bsf_error_val <= digit_10_reach  && digit_10 == 0
            fprintf('The ten %d Gen: %d\n', OFE,gen);%digit_10 = 1;
            outcome(10)=OFE;
            flat_Success=1;
            break;
        end
         
   gen=gen+1;
   
   
   
   %S_bestval 
   PlotActivated=0;
    if (mod(gen,10000)==0 && PlotActivated==1)
        figure(1)
        plot(fitMaxVector(1,:))
        xlabel('Iterations')
        ylabel('Fitness function')
        pause(0.0001) 
    end
    
    Success00_R
  
 
    cont=cont+1;
     
end %---end while ((I_iter < I_itermax) ...
 
%To see when the diversification was used
%  Enhance.std=std_G;
%  Enhance.mean=m_G;
if exist('Act_Evo','var')
  Enhance=Act_Evo;
end
 
if flat_Success==1
    outcome(11)=OFE;
    Fit_and_p=fitMaxVector(1,gen);
else
    outcome(11)=OFE;
    Fit_and_p=fitMaxVector(1,gen-1); %;p2;p3;p4]
end

% subplot(2,1,1)
% plot(mean_Evo)
% subplot(2,1,2)
% semilogy(std_Evo)
 
% VECTORIZED THE CODE INSTEAD OF USING FOR
function pop=genpop(a,b,lowMatrix,upMatrix)
pop=unifrnd(lowMatrix,upMatrix,a,b);

% VECTORIZED THE CODE INSTEAD OF USING FOR
function p=update(p,lowMatrix,upMatrix,BRM,FM_base)
switch BRM
    case 1 %Our method
        %[popsize,dim]=size(p);
        [idx] = find(p<lowMatrix);
        p(idx)=lowMatrix(idx);
        [idx] = find(p>upMatrix);
        p(idx)=upMatrix(idx);
    case 2 %Random reinitialization
        [idx] = [find(p<lowMatrix);find(p>upMatrix)];
        replace=unifrnd(lowMatrix(idx),upMatrix(idx),length(idx),1);
        p(idx)=replace;
    case 3 %Bounce Back
      [idx] = find(p<lowMatrix);
      p(idx)=unifrnd(lowMatrix(idx),FM_base(idx),length(idx),1);
        [idx] = find(p>upMatrix);
      p(idx)=unifrnd(FM_base(idx), upMatrix(idx),length(idx),1);
end

function [FM_ui,FM_base,msg]=generate_trial(method,F_weight, F_CR, FM_pop, FVr_bestmemit,I_NP,I_D,FVr_rot,I_strategyVersion,other)

  FM_popold = FM_pop;                  % save the old population
  FVr_ind = randperm(4);               % index pointer array
  FVr_a1  = randperm(I_NP);                   % shuffle locations of vectors
  FVr_rt  = rem(FVr_rot+FVr_ind(1),I_NP);     % rotate indices by ind(1) positions
  FVr_a2  = FVr_a1(FVr_rt+1);                 % rotate vector locations
  FVr_rt  = rem(FVr_rot+FVr_ind(2),I_NP);
  FVr_a3  = FVr_a2(FVr_rt+1);                
  FM_pm1 = FM_popold(FVr_a1,:);             % shuffled population 1
  FM_pm2 = FM_popold(FVr_a2,:);             % shuffled population 2
  FM_pm3 = FM_popold(FVr_a3,:);             % shuffled population 3
  %FM_mui = rand(I_NP,I_D) < F_CR;  % all random numbers < F_CR are 1, 0 otherwise
  %FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
  
  
  %---- Calculate nearest neighbors for base vectors ----------------
%   for k=1:I_NP
% 	  % Create population with all individuals except k
% 	  FM_popold_without_k = FM_popold([1:k-1,k+1:I_NP],:);
% 	  % Calculate all distances between individual k and the population
% 	  [tmp nnindex]=min(sqrt(sum((ones(I_NP-1,1)*FM_popold(k,:)-FM_popold_without_k).^2,2)));
% 	  % Save the nearest neighbor of individual k
% 	  FM_nn(k,:) = FM_popold(nnindex,:);
%   end
%   % Assign Nearest Neighbors in the base vector
%   FM_pm3 = FM_nn; %For DE/rand1
%   FM_popold=FM_nn;
  %----End: calculate nearest neighbors for base vectors ------------
  
  
  
    if length(F_CR)==1  %Meaning the same F_CR for all individuals
        FM_mui = rand(I_NP,I_D) <= F_CR;  % all random numbers < F_CR are 1, 0 otherwise
        FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
    else %Meaning a different F_CR for each individual
        FM_mui = rand(I_NP,I_D) <= repmat(F_CR,1,I_D);  % all random numbers < F_CR are 1, 0 otherwise
        for i=1:I_NP
            FM_mui (i,randi(I_D))=1;
        end
        FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
    end


    switch method
        case 1,
            FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;     % crossover
            FM_base = FM_pm3;
            msg=' DE/rand/bin';
        case 2,
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            %VEC by FLC
            FM_bm=repmat(FVr_bestmemit,I_NP,1);
            FM_ui = FM_popold + F_weight*(FM_bm-FM_popold) + F_weight*(FM_pm1 - FM_pm2);
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
            FM_base = FM_bm;
            msg=' DE/current-to-best/1';
        case 3, %jDEPerturbated_v3 v4... v7
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            FM_bm=repmat(FVr_bestmemit,I_NP,1);
            if length(F_weight)==1  %Meaning the same F_weight for all individuals
                FM_ui = FM_popold + F_weight*(FM_bm-FM_popold) + F_weight*(FM_pm1 - FM_pm2);
            else
                if  I_strategyVersion==1; %Vortex Algorithm
                        a=other.a;
                        %a = itr / MaxItr; % a value for gammaincinv function
                        %a_par(count)=a;
                        ginv = (1/0.1)*gammaincinv(0.1,a); % compute the new ginv value
                        %ginv_par(count)=ginv;
                        r = ginv * ((other.upperlimit - other.lowerlimit) / 2); %decrease the radius
                        C = r.*randn(I_NP,I_D);
                        FM_ui = bsxfun(@plus, C, FM_bm(1,:));
                end
                
                if  I_strategyVersion==2; %Vortex + HyDE
                     a=other.a; %Linear decrease
                     ginv=exp((1-(1/a^2))); %Exponential decreasing funtion
                     %Not so fast but convergence
                     FM_ui = FM_popold + (repmat(F_weight(:,3),1,I_D).*(FM_pm1 - FM_pm2)) + (ginv)*(repmat(F_weight(:,1),1,I_D).*(FM_bm.*(repmat(F_weight(:,2),1,I_D)+randn(I_NP,I_D))-FM_popold));   % differential variation
                end
                
                if  I_strategyVersion==3; %HyDE (Published in WCCI2018)
                    FM_ui = FM_popold + repmat(F_weight(:,1),1,I_D).*(FM_bm.*(repmat(F_weight(:,2),1,I_D)+randn(I_NP,I_D))-FM_popold) + repmat(F_weight(:,3),1,I_D).*(FM_pm1 - FM_pm2);
                end
                
            end
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
            FM_base = FM_bm;
            msg=' HyDE/current-to-best/1';   
            
    end
return

