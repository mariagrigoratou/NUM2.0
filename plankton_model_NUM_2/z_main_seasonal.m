% Code for the paper "A generic size- and trait-based model of plankton
% communities" 2020 
%
%This code simulates the runs for the SEASONAL ENVIRONEMNT mu_max 
%(i.e. environmental forcing varies with time)
% The mu_max things that need to be defined before runnning is the run-time
% and the number of copepod populations etc (8 between lines 33 to 47)
%if nothing is changed the code will run as it has been run to make the
%figures for the paper
%
% Things to take into account:
%   the model needs to have at least one population of active feeders and one
%   of passive feeders.
%   Most figures have been adapted to the paper, and if the number of state
%   variables is changed some figures might not adapt to it.
%
% Structure of the code:
    % 1- define number of state variables and the environmental forcing
    % 2- the ODE system is solved
    % 3- Some initial figures are created
    % 4- The diagnostics are obtained by re-running the ode function
    % 5- The rest of plots are created
%
% Any questions related to this code can be addressed to me:
% Maria Grigoratou : maria[dot]grigoratou1[at]gmail[dot]com
%
% Maria Grigoratou 26/01/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars

tic
yearsrun=70;      % years of the run

%Parameters setup as in Pompei et al., 2020
nbrP=112;                                 %number of protists size classes
nbrC_act=64;                              %number of copepod active feeders populations
nbrC_pass=24;                             %number of copepods passive feeders populations

P_min=1e-7;                               %minimum size for protists
P_max=0.1;                                %maximum size for protists
Cact_min=0.2;                             %minimum size for active feeders
Cact_max=1000;                            %maximum size for active feeders
Cpass_min=0.2;                            %minimum size for passive feeders
Cpass_max=5;                              %maximum size for passive feeders
C_size_classes=8;                         %Number of size-classes within each copepod population
nbrC=(nbrC_act+nbrC_pass)*C_size_classes; %total number of copepod state variables
nbrF = 1;
%c8 parameters function
param=z_function_parameters(nbrP, P_min, P_max, nbrC_act, Cact_min, Cact_max,...
    nbrC_pass, Cpass_min, Cpass_max, C_size_classes);

% c8 function with predator-prey preferences
[theta_P_P, theta_cop_P, theta_cop_cop, theta_cop]=...
    z_function_feeding_kernels(param);

%c8 function with the seasonal forcing
[T,T_cop,Diff,zmld,light,dzdt]=z_func_seasonal_forcing(param, yearsrun);
Nmax=140; %deep N concentration



%------------------------------------------------
% Other things needed to run the model - Do not modify

% Diagnostics function is not activated
diags=0;
seasonal_switch=1;

global time_counter
time_counter=1;

% indexing for mortality of higher trophic levels      
for i=1:length(param.Wvec)
    idx_bcop{i}= find(param.Wvec>=param.Wvec(i)/10^0.5 & param.Wvec<=param.Wvec(i)*10^0.5);
end
        
tic

%options for the ODE solver
options1 = odeset('Refine',1); 
options2 = odeset(options1,'NonNegative',1:1+param.nbr_P+param.nbr_cops*param.nbr_stages);

%------------------------------------------------
% % Initical conditions as in Pompei 
N0 = zeros(1);
P0 = zeros(1,nbrP);
C0 = zeros(1,nbrC);
 cnt = 1;
  
   for j = 1:cnt
      if j ==1 
        N0 = 1;
        P0(j,:) = 5;
        C0(j,:) =5;

     else
        N0 = N(end);
        P0 = P(end, :);
        C0 = C(end, :);

      end

% Solve ODE
    [t, Y] = ode15s(@z_ode_copepod_model, 1:1:365*yearsrun, [N0, P0, C0], options2, ...
                    param, theta_P_P, theta_cop_P, theta_cop_cop , ...
                    light,T,T_cop,diags,Diff,Nmax,idx_bcop,seasonal_switch,zmld,dzdt);

    toc 

    % Rename output (Y)
    N = Y(:,1);
    P = Y(:,2:1+param.nbr_P);
    C = Y(:,1+param.nbr_P+1:1+param.nbr_P+param.nbr_Ctot);

%% Diagnostics -----------------------------------------------------------
    diags=1;

    Flvl_diag=zeros(size(C));
    dc_diag=zeros(size(C));
    nupos_diag=zeros(size(C));
    pred_C_on_C_diag=zeros(size(C));
    dd_mort_C_diag=zeros(size(C));
    mortP_diag=zeros(size(P));
    pred_P_diag=zeros(size(P));
    pred_C_on_P_diag=zeros(size(P));
    mu_diag=zeros(size(P));
    fg_diag=zeros(size(P));
    fc_diag=zeros(size(C));
    reproduction_diag=zeros(length(t),length(param.Wa));
    fracP_diag=zeros(size(C));
    fracC_diag=zeros(size(C));
    dydt_diag = zeros(size(Y));             % added by MG, March 2021
    gamma_tem_cop_diag = zeros(size(C));    % added by MG, March 2021
    nu_ca_neg_diag = zeros(size(C));        % added by MG, March 2021
    nu_neg_diag = zeros(size(C));           % added by MG, March 2021
    gamma_diag = zeros(size(C));            % added by MG, March 2021
    seeding_C_diag = zeros(size(param.ind_b));        % added by MG, March 2021
    MLD_diag = zeros(size(t));
    MLD_C_diag = zeros(size(C)); 
    I_env_diag = zeros(size(C));
    k_diag = zeros(size(C));

    for i=1:length(t)

          y_init=[N(i), P(i,:), C(i,:)]';
          [dydt, Flvl, d_c, pred_C_on_C, dd_mort_C, nu_pos, mortP , pred_P, pred_C_on_P,mu,fg,fc,...
              reproduction,fracP,fracC, gamma_temp_cop,nu_ca_neg, nu_neg, gamma, seeding_C, MLD, MLD_C, I_env, k] = z_ode_copepod_model(i, y_init,...
                  param, theta_P_P, theta_cop_P, theta_cop_cop , ...
                   light,T,T_cop, diags,Diff,Nmax,idx_bcop,seasonal_switch,zmld,dzdt); 


          Flvl_diag(i,:)=Flvl;     
          dc_diag(i,:)=d_c;   
          nupos_diag(i,:)=nu_pos;
          pred_C_on_C_diag(i,:)=pred_C_on_C;
          dd_mort_C_diag(i,:)=dd_mort_C;
          mortP_diag(i,:)=mortP;
          pred_P_diag(i,:)=pred_P;
          pred_C_on_P_diag(i,:)=pred_C_on_P;
          mu_diag(i,:)=mu;
          fg_diag(i,:)=fg;
          fc_diag(i,:)=fc;
          reproduction_diag(i,:)=reproduction;
          fracP_diag(i,:)=fracP; 
          fracC_diag(i,:)=fracC;
          dydt_diag(i,:) = dydt;                             % added by MG, March 2021
          gamma_tem_cop_diag(i, :) = gamma_temp_cop;         % added by MG, March 2021
          nu_ca_neg_diag(i,:) = nu_ca_neg;                   % added by MG, March 2021
          nu_neg_diag(i,:) = nu_neg;                         % added by MG, March 2021
          gamma_diag(i,:) = gamma;                           % added by MG, March 2021
          seeding_C_diag(i,:) = seeding_C;                   % added by MG, March 2021
          MLD_diag(i) = MLD;                                 % added by MG, March 2021
          MLD_C_diag(i,:) = MLD_C;                           % added by MG, March 2021
          I_env_diag(i, :) = I_env;                          % added by MG, March 2021
          k_diag(i, :) = k;                                  % added by MG, March 2021

    end

    Flvl_mean=mean(Flvl_diag(ceil(t(end)*3/4):end,:),1);
    Flvl_max=max(Flvl_diag(ceil(t(end)*3/4):end,:));
    Flvl_min=min(Flvl_diag(ceil(t(end)*3/4):end,:));

    nu_mean=mean(nupos_diag(ceil(t(end)*3/4):end,:),1);
    dc_mean=mean(dc_diag(ceil(t(end)*3/4):end,:),1);
    predCC_mean=mean(pred_C_on_C_diag(ceil(t(end)*3/4):end,:),1);
    dd_mort_C_mean=mean(dd_mort_C_diag(ceil(t(end)*3/4):end,:),1);
    mortP_mean=mean(mortP_diag(ceil(t(end)*3/4):end,:),1);
    pred_P_mean=mean(pred_P_diag(ceil(t(end)*3/4):end,:),1);
    pred_C_on_P_mean=mean(pred_C_on_P_diag(ceil(t(end)*3/4):end,:),1);
    mu_mean=mean(mu_diag(ceil(t(end)*3/4):end,:),1);
    fg_mean=mean(fg_diag(ceil(t(end)*3/4):end,:),1);
    fc_mean=mean(fc_diag(ceil(t(end)*63/4):end,:),1);

    Pmean=mean(P(ceil(t(end)*3/4):end,:),1);
    Cmean=mean(C(ceil(t(end)*3/4):end,:),1);
    deltaC=param.deltaC;
    
 end
%% figures
% % %
figure
subplot(2,3,1)
plot(t,N,'k')

subplot(2,3,2)
plot(t,P)
legend
title('protists')

subplot(2,3,3)
plot(t,C(:,param.ind_a(1:nbrC_act)))
legend
title('active cop')

subplot(2,3,4)
plot(t,C(:,param.ind_a(nbrC_act+1:end)))
legend
title('passive cop')
%
figure
subplot(2,3,1)
plot(t(end-365:end),N(end-365:end),'k')

subplot(2,3,2)
plot(t(end-365:end),P(end-365:end,:))
title('protists')

subplot(2,3,3)
plot(t(end-365:end),C(end-365:end,param.ind_a(1:nbrC_act)))
%legend
title('active cop')

subplot(2,3,4)
plot(t(end-365:end),C(end-365:end,param.ind_a(nbrC_act+1:end)))
%legend
title('passive cop')

% %Main plot seasonal environment

yr=1;
yrend=0*365; 


fsize=10;
rownbs=1;

s1=1e-6;
s2=1e-5;
s3=1e-4;
s4=1e-3;
s5=1e-2;

tinterval=t(end-364*yr):30.5:t(end)-yrend; %in Pompei et al is end - 365. I changed to 364 to get the last 365 days - MG Mar 2021
Gyear=P(end-729*yr:end-yrend,:);  %in Pompei et al is end - 365. I changed to 364 to get the last 365 days - MG Jan 2021
gg1=sum(Gyear(:,param.V<s1),2);
gg2=sum(Gyear(:,param.V>=s1 & param.V<s2),2);
gg3=sum(Gyear(:,param.V>=s2 & param.V<s3),2);
gg4=sum(Gyear(:,param.V>=s3 & param.V<s4),2);
gg5=sum(Gyear(:,param.V>=s4 & param.V<s5),2);
gg6=sum(Gyear(:,param.V>=s5),2);

gg=cat(2,gg1,gg2,gg3,gg4,gg5,gg6);
months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
bluemap=(brewermap(6,'YlGnBu'));


x0=0;
y0=0;
width=17;
height=16;
grey=[0.5, 0.5, 0.5];
yellow=[244, 185, 66]./255;

cp= [0.5882    0.5882    0.5882];

 %%
% %size classes to make body-mass ranges of copepods
sc1=1e-2;
sc2=1e-1;
sc3=1e0;
sc4=1e1;
sc5=1e2;

% %group copepods in the body-mass ranges
Cyear=C(end-729*yr:end-yrend,param.ind_act); % in Pompei et al 2021 was end - 365. I changed it to 364 to get the last 365 days.- MG Jan 2021
cc1_a=sum(Cyear(:,param.W(param.ind_act)<sc1),2);
cc2_a=sum(Cyear(:,param.W(param.ind_act)>=sc1 & param.W(param.ind_act)<sc2),2);
cc3_a=sum(Cyear(:,param.W(param.ind_act)>=sc2 & param.W(param.ind_act)<sc3),2);
cc4_a=sum(Cyear(:,param.W(param.ind_act)>=sc3 & param.W(param.ind_act)<sc4),2);
cc5_a=sum(Cyear(:,param.W(param.ind_act)>=sc4 & param.W(param.ind_act)<sc5),2);
cc6_a=sum(Cyear(:,param.W(param.ind_act)>=sc5),2);

Cyear_p=C(end-729*yr:end-yrend,param.ind_pass); % in Pompei et al 2021 was end - 365. I changed it to 364 to get the last 365 days.- MG Jan 2021
cc1_p=sum(Cyear_p(:,param.W(param.ind_pass)<sc1),2);
cc2_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc1 & param.W(param.ind_pass)<sc2),2);
cc3_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc2 & param.W(param.ind_pass)<sc3),2);
cc4_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc3 & param.W(param.ind_pass)<sc4),2);
cc5_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc4 & param.W(param.ind_pass)<sc5),2);
cc6_p=sum(Cyear_p(:,param.W(param.ind_pass)>=sc5),2);

cc_a=cat(2,cc1_a,cc2_a,cc3_a,cc4_a,cc5_a,cc6_a);
cc_p=cat(2,cc1_p,cc2_p,cc3_p,cc4_p,cc5_p,cc6_p);

% % %save the values in a mat file -- % added by MG, March 2021
% save(['C:\Users\mgrigoratou\Documents\GitHub\Copepod_sizebased_model-master_no_fp\plankton_model_v2\final_march22\model_output\\model_output_mld_mab_slope_2000_2010_CAMILA_SETUP_summer' num2str(j) '.mat'],'C', 'P', 'N', 'C0', 'P0', 'N0','fc_mean', 'fg_mean',...
%         'nu_mean', 'pred_C_on_P', 'pred_P_mean', 'mortP_mean','dd_mort_C_mean',...
%         'predCC_mean', 'dc_diag', 'nu_mean','dc_mean', 'nu_mean', 'Flvl_diag', 'dc_diag', 'nupos_diag', ...
%         'pred_C_on_C_diag', 'dd_mort_C_diag', 'mortP_diag','pred_P_diag', 'pred_C_on_P_diag', ...
%          'mu_diag','fg_diag', 'fc_diag', 'reproduction_diag','fracP_diag','fracC_diag','dydt_diag', ...
%          'gamma_tem_cop_diag','nu_ca_neg_diag','nu_neg_diag','gamma_diag','seeding_C_diag',...
%          'MLD_diag','MLD_C_diag','I_env_diag','k_diag', 'param', 'nbrC_act', 'cc_a', 'cc_p', 'gg');
%      
% 
% % %save only the biomass output for Heatwave scenarios -- % added by MG, March 2021
% C_4dg = C; 
% P_4dg = P;
% N_4dg = N; 
% 
% cc_a_warming = cat(2,cc1_a,cc2_a,cc3_a,cc4_a,cc5_a,cc6_a);
% cc_p_warming = cat(2,cc1_p,cc2_p,cc3_p,cc4_p,cc5_p,cc6_p);
% gg_warming = cat(2,gg1,gg2,gg3,gg4,gg5,gg6);
% 
% save(['C:\Users\mgrigoratou\Documents\GitHub\Copepod_sizebased_model-master_no_fp\plankton_model_v2\final_march22\model_output\\model_output_mld_mab_slope_2000_2010_CAMILA_SETUP_summer' num2str(j) '.mat'],...
%         'C', 'N', 'P');
%      
     

