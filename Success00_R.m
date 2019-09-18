%Lock like this
%Optimize Happy cat function

%   deParameters.I_strategy=3;
%     deParameters.I_strategyVersion=2;
%     deParameters.I_itermax_DE=1e5; %Duration of vortex function
%     
%     deParameters.I_itermax= 5e6; %Maximum budget
% 
%     deParameters.I_NP=50; 100 with happy cat
%     deParameters.F_weight=0.5;
%     deParameters.F_CR=0.5;
%     deParameters.I_bnd_constr=3;

%Funct 7: 
%Funct 8: 
%Funct 9: 2 digits


%Enhanced population with a simple mechanism
if ~exist('countStagnation','var'),countStagnation=0;end

% if ~exist('m_G','var'),m_G=zeros(I_itermax,I_D);end
% if ~exist('std_G','var'),std_G=zeros(I_itermax,I_D);end
 if ~exist('act_cont','var'),act_cont=1;end
 if ~exist('Act_Evo','var'),Act_Evo=0;end
% 
% std_G(gen,:)=std(FM_pop);
% m_G(gen,:)=mean(FM_pop);
%Act_Evo(gen)=0;

if sum(ind)==0 %FailedIteration
    countStagnation = countStagnation + 1;
else
    countStagnation = 0;
end


if countStagnation==1e4
     Act_Evo(act_cont)=gen;
     act_cont=act_cont+1;
        % do_smt=1;
        countStagnation=0;

        [~,ind]=sort( S_val);
      usa=10;
      cont1=1;
      for ii=1:usa
%      FM_bestC=FM_pop(:,ind(1:5)); %Select the best percentage
      micro=(FM_pop(ind(ii),:)-FVr_minbound)./(FVr_maxbound-FVr_minbound);
      sigma=10e-4;
      randnormal=normrnd(repmat(micro,I_NP/usa,1),repmat(sigma,I_NP/usa,I_D));
       
      FM_pop(cont1:cont1-1+I_NP/usa,:)= minPositionsMatrix(1:I_NP/usa,:)+(maxPositionsMatrix(1:I_NP/usa,:)-minPositionsMatrix(1:I_NP/usa,:)).*randnormal;
      cont1=cont1+I_NP/usa;
      end
      
        f = cec19_func(FM_pop',fnc);
        S_val=f';
        OFE=OFE+I_NP;%+I_NP;
        %Retrieve best individual information
        S_val(I_best_index) = S_bestval;
        FM_pop(I_best_index,:)=FVr_bestmemit; %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      
        
       
        
%           F_weight_old=repmat(deParameters.F_weight,I_NP,3);
%           F_weight= F_weight_old;
%           F_CR_old=repmat(deParameters.F_CR,I_NP,1);
%           F_CR=F_CR_old;

end


if mod(gen,I_itermax_DE)==0
   cont=1; %Reset the DE_vortex function
end

 
