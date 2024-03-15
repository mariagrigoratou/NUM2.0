function [dydt, F_lvl, d_c, pred_C_on_C, dd_mort_C, nu_pos, mortP, pred_P, pred_C_on_P,mu,fg,fc,...
    reproduction,fracP,fracC, gamma_temp_cop, nu_ca_neg, nu_neg, gamma, seeding_C, MLD, MLD_C, I_env, k]=...
    z_ode_copepod_model(t, y,...
    param,theta_P, theta_cop_P, theta_cop_cop, Ls,T_inn,T_inn_cop,diags,Diff, Nmax, idx_bcop,...
    seasonal_switch,zmld2,dzdt2)
% Code for the paper "A generic size- and trait-based model of plankton
% communities" 2020 
%
%This code has all the ODE of the model
% parameters and variables have different names from the ones in the 
% paper. The code is set up with the control enviromental conditons. Go to 
% z_func_seasonal_forcing to set up the Heatwave scenario you want to test.  

% The INPUTS to this function are:
    %param: parameters
    %theta_P: preference of protists for protists
    %theta_cop_P: preference of copepods for protists
    %theta_cop_cop: preference of copepods for copepods
    %theta_cop_F: preference of copepods for fecal pellets
    %Ls: surface light,
    %T_inn: temperature
    %GRIGORATOU TEMPERATURE
    % t_factor           : normalization factor for temp, unitless
    % a_temp             : temp coefficient , Kelvin
    % t_ref              : temp reference , Kelvin
    % b_width            : factor determining the width of norms, 1/K
    % b                  : decay coefficient, unitless 
    % t_opt              : norm optimumn temp, Kelvin
    % temp               : ambient temp
    % gamma_temp         : temp parameter in growth
    % T_inn_kelvin       : ambient temperature in Kelvin (calls the temp parameter from z_func_seasonal_forcing) 
    % T_inn_kelvin_hw    : ambient temperature in Kelvin for the
    % experiments when the direct effects of temp are on copepods or
    % protists only. calls the temp_cop parameter from z_func_seasonal_forcing) 
    
    %diags: defines whether we want to have the diagnostics as output. 
        %if diags=0 the ODE system is solved, if diags=1 we use the function to
        %get the diagnostics
    %Diff: input of nitrogen in the system (rho in the paper)
    %Nmax: deep nutrients
    %idx_bcop: indexes for the mortality by higher trophic levels (saves time to define them outside)
    %seasonal_switch: switches between constant and seasonal environment ()
        %seasonal_switch=0 --> constant environment
        %seasonal_switch=1 --> seasonal environment
    %zmld2: mixed layer depth
    %dzdt2: rate of change of the mixed layer
        
% Structure of the code:
    % 1- defines environemntal inputs according to the constant/seasonal
    % scenario
    % 2- system of ODES
    % 3- Groups ODE outputs (dydt) or diagnostics
    % 4- Checks for mass balance
%
%
% Any questions related to this code can be addressed to me:
% cam.serra90@gmail.com
%
% Camila Serra-Pompei 25/08/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<14
    seasonal_switch=0;    
end

y=y';

global time_counter

%checks of there are NaNs
if any(isnan(y))
      disp('there are NaNs inside the ODE!!=');
end  
%checks of there are imaginary numbers
if any(imag(y))
      disp('there are NaNs inside the ODE!!=');
end  
%speeds-up runnin time and small values going below 0
y=max(y,1e-150);


%set state variables
N = y(1);
P = y(2:1+param.nbr_P);
C = y(1+ param.nbr_P+1:1+ param.nbr_P + param.nbr_Ctot);
C=max(C,0);


% define environmental parameters for each scenario
if seasonal_switch==1                                  % if seasonal environment
    if t==0
        Is=Ls(1);                                      % Light at surface
        dzdt=dzdt2(1);                                 % Rate of change of ML
        zmld=zmld2(1);                                 % MLD
        T_inn=T_inn(1);                                % Temperature
        T_inn_cop = T_inn_cop(1);                      % temp for cop 
        Diff=((param.Diff_min + max(0,dzdt))./(zmld)); % Input of nitrogen
        
    else
        Is=Ls(ceil(t));
        dzdt=dzdt2(ceil(t));
        zmld=zmld2(ceil(t));
        T_inn=T_inn(ceil(t));
        T_inn_cop=T_inn_cop(ceil(t));
        Diff=((param.Diff_min + max(0,dzdt))./(zmld));

    end
    
    ktot = param.kw + param.kchl .* sum(P);
    L = (Is./(ktot*zmld)).*(1-exp(-ktot.*zmld));        % Mean PAR in surface layer
    MLD_P= - ((param.Diff_min + max(0,dzdt))./zmld).*P; %effects of ML dynamics on protists
    MLD = -dzdt/zmld;
    MLD_C= - (dzdt/zmld) .* C;                          %effects of ML dynamics on copepods
    seeding_P = param.flow'.*(param.inputP - P);
    seeding_C = param.flowC'.*(param.inputC - C(param.ind_b));
    
else % if constant environment
    
  MLD_P=0;
  MLD_C=0;
  seeding_P=0;
  seeding_C=0;
  L=Ls;
  
end

temp_option = 'Grigoratou'; % 'Pompei' (Q10) or 'Grigoratou' (MIT- Darwin) temperature 

switch temp_option 
    case 'Pompei'                                                  % temp units: oC         
        I=param.I'.*param.Q_I.^((T_inn-param.Tref)/10);            % IR cops
        alpha_N=param.alpha_N.*param.Q_N.^((T_inn-param.Tref)/10); % protist
        mu_max=param.mu_max.*param.Q_k.^((T_inn-18)/10);           % protist
        remin=param.remin.*param.Q_k.^((T_inn-param.Tref)/10);     % environment
        alpha_F=param.alpha_F.*1.5.^((T_inn-param.Tref)/10);       % protist
        alpha_c=param.alpha_c'.*1.5.^((T_inn-param.Tref)/10);      % CR active & passive
        VN=param.VN.*2.^((T_inn-18)/10);                           % protist
        VF=param.VF.*2.^((T_inn-18)/10);                           % protist
        VL=param.VL.*2.^((T_inn-18)/10);                           % protist
        respP=param.R.*param.Q_I.^((T_inn-param.Tref)/10);         % protist
               

         
    case 'Grigoratou'                         %added by MG Dec 2020, temp units: Kelvin
        t_factor = 0.8;                       % normalization factor for temp function
        a_temp = -4000;                       % 
        t_ref = 288.15;                       % temperature reference 288.15 (15oC) Camila's 
        b_width = 3*10^(-4);                  % width of norms, default: 3*10^(-4) 
        b_temp = 4;                           % decay coefficient for norms
        T_inn_kelvin = T_inn +273.15;         % ambient temp in Kelvin --> calls the temp from the z_func_seasonal_forcing this is the set up for HW on both protists and copepods
        T_inn_kelvin_HW = T_inn_cop + 273.15; % temp for testing the HWs. Added in both prot and cops or only at one of them depending on the tested hypothesis
        % temp limitation for protists - added by MG Dec2020
        gamma_temp_P = t_factor * exp(a_temp * (1/T_inn_kelvin - 1/t_ref))...
            *exp(-b_width * abs(T_inn_kelvin - param.temp_opt_P).^b_temp);
              
        % temp limitation for copepods - added by MG Dec2020
        gamma_temp_cop = t_factor * exp(a_temp * (1/T_inn_kelvin - 1/t_ref))...
            *exp(-b_width * abs(T_inn_kelvin - param.temp_opt_cop).^b_temp);

        
        remin=param.remin;   %reminalization 
        % %protists
        alpha_N=param.alpha_N.*gamma_temp_P;       
        mu_max=param.mu_max.* gamma_temp_P;        
        VN=param.VN.* gamma_temp_P;             
        VF=param.VF.* gamma_temp_P;             
        VL=param.VL.* gamma_temp_P;             
        alpha_F=param.alpha_F.* gamma_temp_P; 
        respP=param.R.* max(gamma_temp_P); %%% add max like I_env for copepods, remove max if you don't want it like I_env for cops
        % %cops
        alpha_c=param.alpha_c'.* gamma_temp_cop'; 
        I = param.I'.* gamma_temp_cop'; 
        I_env = param.I'.*max(gamma_temp_cop);
        
        
    otherwise
    disp ('wrong option')
    dbstop


end
%transpose to the right dimensions (small corrections)
alpha_c=alpha_c';


% % protists---------------------------------------------------------------

%food for protists
E_P_tot=P*theta_P';
%food for Copepods
E_C_tot=(theta_cop_P*P')' + (theta_cop_cop*C')';

%Uptake of resources for protists
J_F=(VF.*alpha_F.*E_P_tot)./(alpha_F.*E_P_tot+VF);
J_N=(VN.*param.Qcn.*alpha_N.*N)./(param.Qcn.*alpha_N.*N+VN);
J_L=(VL.*param.alpha_L.*L)./(param.alpha_L.*L+VL);

%respiration
J_R=respP;

%growth rate
mu=min(J_L + J_F - J_R , J_N+ J_F);

%Nleaks
leaks=max(0,J_N - J_L + J_R);

%feeding level protists
fg=(alpha_F.*E_P_tot)./(alpha_F.*E_P_tot+VF);

%predation by protists on protists
pred_P=(J_F.*P./E_P_tot)*theta_P;

% %copepods ---------------------------------------------------------------
%Feeding level of copepods
F_lvl=(alpha_c'.*E_C_tot)./(alpha_c'.*E_C_tot+I).*gamma_temp_cop';
k=param.eff.*I_env'.*(0.1+0.2.*F_lvl');
fc=k'./(param.eff'.*I); %critical feeding level

%predation by copepods on protists
pred_C_on_P=(I.*F_lvl.*C./E_C_tot)*theta_cop_P;

%background mortality of protists
mortP=mu_max.*0.03.*P./param.ratio_V;

% Protists ODE
dPdt = mu.*P .*gamma_temp_P-  mortP.*P - pred_P.*P - pred_C_on_P.*P + seeding_P + MLD_P;

% %background mortality copepods
mmax=(param.no_small'.*I_env.*0.003)./param.deltaratio(:)'; 
dd_mort_C=zeros(size(C));
for i=1:length(param.Wvec)
   Bcop=sum(C(:,idx_bcop{i}));  
   if imag(Bcop)
       
   end
   dd_mort_C(:,i)=mmax(i).*C(i).^(0.2).*Bcop.^(1-0.2);   
   if any(imag(dd_mort_C))
       
   end
end
% 

 
% Net energy gain for all stages
nu=param.migration'.*param.eff'.*I.*(F_lvl-fc); 
nu_pos=max(0,nu);       % positive growth
nu_neg=-min(0,nu);      % starvation mortality
nu_ca_neg= min(0,nu);   % starvation of the adults

%Predation from copepods on copepods
pred_C_on_C=(I.*F_lvl.*C./E_C_tot)*theta_cop_cop;
    

%maturation rates
d_c=pred_C_on_C+nu_neg+dd_mort_C; % total mortality
gamma=zeros(size(C));             % describes the transfer of biomass between size classes, i.e.somatic growth
gamma(param.ind_j)=(nu_pos(param.ind_j)-d_c(param.ind_j))./(1-param.z(param.ind_j)'.^(1-(d_c(param.ind_j)./nu_pos(param.ind_j))));

    

%ODE copepod
dCdt=zeros(size(C));

% Size at birth

dCdt(param.ind_b) = param.rep_eff.* nu_pos(param.ind_a).*C(param.ind_a) ...
                    +(nu(param.ind_b)...
                     - gamma(param.ind_b)...
                     - pred_C_on_C(param.ind_b) - dd_mort_C(param.ind_b)).*C(param.ind_b) + seeding_C ;
                 
                          
% In-between size-classes
dCdt(param.ind_rest)= gamma(param.ind_rest-1).*C(param.ind_rest-1) ...
                       + (nu(param.ind_rest)...
                       - gamma(param.ind_rest)...
                       - pred_C_on_C(param.ind_rest) - dd_mort_C(param.ind_rest)).*C(param.ind_rest) ;

% Adults  - as in Pompei                               
dCdt(param.ind_a) = gamma(param.ind_a-1).*C(param.ind_a-1)...
                     + (nu_ca_neg(param.ind_a) ...
                     - pred_C_on_C(param.ind_a) - dd_mort_C(param.ind_a)).*C(param.ind_a);

                
dCdt= dCdt+ MLD_C; %effects from ML dynamics  
 

% all dead organisms (rate) fraction to be remineralised
stuff_dying=sum(mortP.*P) + sum(dd_mort_C.*C);  

%Nitrogen ode
dNdt=Diff.*(Nmax-N) + ( - sum(J_N .* P, 2) + sum((leaks.*P),2) +  sum(k'.*C,2)...
         + sum((1-param.rep_eff).* nu_pos(param.ind_a).*C(param.ind_a),2) ...
         + param.remin_frac*(stuff_dying))./param.Qcn;
             

% display time
if seasonal_switch == 1
       
    %time counter (displays time)
    if mod(ceil(t),365)==0 
        if time_counter~=mod(ceil(t),1000)
            disp(['Year: ' num2str(ceil(t)/365)]) 
            time_counter=mod(ceil(t),1000);
        end
        
    else
        time_counter=1; 
    end
    
    
    if diags==1 %some more diagnostics
       reproduction=param.rep_eff.* nu_pos(param.ind_a);
       fracP=(theta_cop_P*P')'./E_C_tot;
       fracC=(theta_cop_cop*C')'./E_C_tot;


    end
    
else
    
    if diags==1 %some more diagnostics
       reproduction=param.rep_eff.* nu_pos(param.ind_a);
       fracP=(theta_cop_P*P')'./E_C_tot;
       fracC=(theta_cop_cop*C')'./E_C_tot;
    end
             
    %time counter (displays time)
    if mod(ceil(t),1000)==0 
        if time_counter~=mod(ceil(t),1000)
            disp(['Days: ' num2str(ceil(t))]) 
            time_counter=mod(ceil(t),1000);
        end
        
    else
        time_counter=1; 
    end

end
            

% Check for mass balance ------------------------------
if t<10000 %do it only in first time-steps to save time   
    if seasonal_switch==1
        dNtotdt=dNdt + (sum(dPdt) + sum(dCdt) )./param.Qcn;
        dNindt=Diff.*Nmax + (sum(seeding_P) + sum(seeding_C))./param.Qcn;
        dNoutdt=Diff.*N + (1/param.Qcn)*((1-param.remin_frac)*(stuff_dying) ...
                - sum(MLD_P) - sum(MLD_C));
         massbalance=dNtotdt-dNindt+dNoutdt;
                
    else
        dNtotdt=dNdt + (sum(dPdt) + sum(dCdt))./param.Qcn;
        dNindt=Diff.*Nmax;
        dNoutdt=Diff.*N + (1/param.Qcn)*((1-param.remin_frac)*(stuff_dying) ...
                + sum((param.sink_fp./10).*F));
        massbalance=dNtotdt-dNindt+dNoutdt;

    end
    
    if massbalance>1e-10
        disp('No mass balance!') 
        disp ('massbalance')
        massbalance 
        
    end

end
 

dydt = [dNdt; dPdt'; dCdt'];

if any(isnan(dydt))
      disp('there are NaNs inside the ODE!!=');
end  

if any(imag(dydt))
      disp('there are NaNs inside the ODE!!=');
end
% save all the Ode steps -added by MG March 2021 
% is very slow, run it only for 1 year if you want to save the values you are testing. 
% if exist('temp.mat', 'file')
%      next_gamma_temp_cop = gamma_temp_cop;
%      load temp.mat gamma_temp_cop
%      gamma_temp_cop = cat(2, next_gamma_temp_cop, gamma_temp_cop);
%      save temp.mat gamma_temp_cop -append
%      gamma_temp_cop = next_gamma_temp_cop;
% else
%     save temp.mat gamma_temp_cop
% end

% if exist('dc.mat', 'file')
%      next_dCdt = dCdt;
%      load dc.mat dCdt
%      dCdt = cat(2, next_dCdt.', dCdt);
%      save dc.mat dCdt -append
%      dCdt = next_dCdt;
% else
%     save_dCdt = dCdt.';
%     save dc.mat save_dCdt
% end
% test end - added by MG March 2021

end
 
