function param=z_function_parameters(P_size_classes, P_min, P_max, Cact_populations, Cact_min, Cact_max,...
    Cpass_populations, Cpass_min, Cpass_max, C_size_classes, nbrD,PPMR)

% Here we create the whole set of parameters needed to run the model
%note that te names of parameters here are different than in the paper or
%the ones you introduced initially. But everything has been commented. If
%you have any questions feel free to email me: mcsp@aqua.dtu.dk

%Camila Serra 12/03/2020
%-------------------------------------------------------------------

% %Make grid of protists
% param.nbr_P=P_size_classes; %number of size-classes
% V_points=param.nbr_P+1; %number of points in the grid
% param.min_P=log10(P_min); %minimum size generalists in log10
% param.max_P=log10(P_max); %maximum size generalists in log10
% param.V_grid=logspace(param.min_P,param.max_P,V_points); %generalists grid
% param.V_up=param.V_grid(2:end); %upper sizes
% param.V_dw=param.V_grid(1:end-1); %lower sizes
% param.V=geomean([param.V_up;param.V_dw]); %central sizes
% param.delta_V=(param.V_up)-(param.V_dw); %bin size
% param.ratio_V=(param.V_up./param.V_dw); %bin size
% 
% %Make grid of copepods
% C_sp_act=Cact_populations; %number of species (populations) active feeders
% C_sp_pass=Cpass_populations;%5; %number of species (population) passive feeders
% param.C_sp_act=C_sp_act;
% param.C_sp_pass=C_sp_pass;
% param.C_sp=C_sp_act+C_sp_pass; %total number of populations
% param.nbr_discr=C_size_classes;
% param.min_cop_act=log10(Cact_min); % minimum size of adult active copepods in log10
% param.max_cop_act=log10(Cact_max); % maximum size of adult active copepods in log10
% param.min_cop_pass=log10(Cpass_min); % minimum size of adult passive copepods in log10
% param.max_cop_pass=log10(Cpass_max); % maximum size of adult passive copepods in log10
% C_grid_act=zeros(param.nbr_discr,C_sp_act); % each column is a population of active copepods, 
%                                             % each row is a size bin, from juveniles (top) to adults (bottom)
% C_grid_pass=zeros(param.nbr_discr,C_sp_pass); %same as above for passive copepods
% Ca_sp_act=logspace(param.min_cop_act,param.max_cop_act,C_sp_act); %vector of adult copepod size
% Ca_sp_pass=logspace(param.min_cop_pass,param.max_cop_pass,C_sp_pass); %vector of adult copepod size
% 
% 
size_option = 'Grigoratou';
switch size_option 
    case 'Grigoratou'
        M = readtable('C:\Users\mgrigoratou\Documents\GitHub\Copepod_sizebased_model-master_no_fp\plankton_model_v2\\species_dataframe_8_temp.xlsx');
        % %extract arrays from the table
        param.pft = (table2array(M(:,{'pft_numb'})))';
        param.nbr_total = height(M(:,{'pft_numb'}));
        param.min_size = (table2array(M(:, {'min_size'})))'; % transfer vertical to horizontal array
        param.max_size = (table2array(M(:, {'max_size'})))';
        param.temp_opt = (table2array(M(:, {'temp_opt'})))';

        % %define zero arrays for organisms
        %protists
        param.size_bin_P = 1;              % protists' size bin
        param.nbr_P = 0;                   % number of protists in the dataframe
        param.min_P = zeros();             % min size of protists in log10
        param.max_P = zeros();             % max size of protists in log10
        param.temp_opt_P = zeros();        % temp optima of protists 
        %cops 
        param.nbr_discr = 8;               % C_size_classes 
        %active cops
        C_sp_act = 0;                      % number of active cops in the dataframe
        param.min_cop_act = zeros();       % min size of active cops in log10
        param.max_cop_act = zeros();       % max size of active cops in log10
        param.temp_opt_cop_act = zeros();  % temp optima of active cops
        %passive cops
        C_sp_pass = 0;                     % number of passive cops in the dataframe
        param.min_cop_pass = zeros();      % min size of passive cops in log10
        param.max_cop_pass = zeros();      % max size of passsive cops in log10
        param.temp_opt_cop_pass = zeros(); % temp optima of passive cops


        for z= 1:length(param.pft)
            switch param.pft(z)
             case 1
                %disp('protist');
                param.nbr_P = param.nbr_P +1; %number of size-classes
                param.min_P(z) = param.min_size(z);
                param.max_P(z) = param.max_size(z);
                param.temp_opt_P(z) = param.temp_opt(z);
            case 2
                %disp('active copepod');
                C_sp_act = C_sp_act + 1;
                param.min_cop_act(z) = param.min_size(z);
                param.max_cop_act(z) = param.max_size(z);     
                param.temp_opt_cop_act(z) = param.temp_opt(z); 
            case 3
                %disp('passive copepod');
                C_sp_pass = C_sp_pass + 1;
                param.min_cop_pass(z) = param.min_size(z);
                param.max_cop_pass(z) = param.max_size(z);     
                param.temp_opt_cop_pass(z) = param.temp_opt(z);
          
                otherwise
                    disp('error: organism unknown')
                    dbstop
            end
        end

        % %remove zeros from the arrays
        % %protists
        param.min_P = nonzeros(param.min_P)';
        param.max_P = nonzeros(param.max_P)'; 
        param.temp_opt_P = nonzeros(param.temp_opt_P)'; 
        % %active copepods
        param.min_cop_act = nonzeros(param.min_cop_act);
        param.max_cop_act = nonzeros(param.max_cop_act);
        param.temp_opt_cop_act = nonzeros(param.temp_opt_cop_act);
        Ca_sp_act = param.max_cop_act.'; % all max sizes - adult copepods
        % %passive copepods
        param.min_cop_pass = nonzeros(param.min_cop_pass);
        param.max_cop_pass = nonzeros(param.max_cop_pass);
        param.temp_opt_cop_pass = nonzeros(param.temp_opt_cop_pass);
        Ca_sp_pass = param.max_cop_pass.'; % all max sizes - adult copepods
        % 
        % %temperature for all size bins for cops
        param.temp_opt_cop_act = repmat(param.temp_opt_cop_act.',param.nbr_discr, 1);
        param.temp_opt_cop_pass = repmat(param.temp_opt_cop_pass.',param.nbr_discr, 1);
        param.temp_opt_cop = [param.temp_opt_cop_act,param.temp_opt_cop_pass];
        param.temp_opt_cop = param.temp_opt_cop(:);

        %
        % %grid of protists
        V_points = param.nbr_P + 1;
        nbrP = param.nbr_P;
        param.V_grid = [param.min_P(1:14), param.max_P(14)]; 
        param.V = geomean([param.min_P;param.max_P])    % central sizes
        param.V_up = param.V_grid(2:15);                % upper sizes
        param.V_dw = param.V_grid(1:15-1);              % lower sizes
        param.delta_V = (param.V_up)-(param.V_dw);      % bin size
        param.ratio_V = (param.V_up./param.V_dw);       % bin size
        param.prot_set = (nbrP/14);                     % set of protists (14 groups) per dataframe - use that for param.ratio_V and delta_V
        param.ratio_V = repmat(param.ratio_V.',param.prot_set,1.').';
        param.delta_V  = repmat(param.delta_V.', param.prot_set, 1).';
        param.V_up = repmat(param.V_up.',param.prot_set,1.').';        % upper sizes
        param.V_dw = repmat(param.V_dw.',param.prot_set,1.').';        % lower sizes

        % %cops total number of populations
        C_grid_act = zeros(param.nbr_discr,C_sp_act);
        C_grid_pass = zeros(param.nbr_discr,C_sp_pass);  
        C_sp = C_sp_act + C_sp_pass;                % total number of populations
        for i=1:C_sp_act
            C_grid_act(:,i)=logspace(log10(param.min_cop_act(i)),(log10(param.max_cop_act(i))),param.nbr_discr); %copepod grid 
        end
        for i=1:C_sp_pass
            C_grid_pass(:,i)=logspace(log10(param.min_cop_pass(i)),(log10(param.max_cop_pass(i))),param.nbr_discr); %copepod grid 
        end
                   
    case 'Pompei'        
        %Make grid of protists
        param.nbr_P=P_size_classes; %number of size-classes
        V_points=param.nbr_P+1; %number of points in the grid
        param.min_P=log10(P_min); %minimum size generalists in log10
        param.max_P=log10(P_max); %maximum size generalists in log10
        param.V_grid=logspace(param.min_P, param.max_P, V_points); %generalists grid
        param.V_up=param.V_grid(2:end); %upper sizes
        param.V_dw=param.V_grid(1:end-1); %lower sizes
        param.V=geomean([param.V_up;param.V_dw]); %central sizes
        param.delta_V=(param.V_up)-(param.V_dw); %bin size
        param.ratio_V=(param.V_up./param.V_dw); %bin size
        %Make grid of copepods
        C_sp_act=Cact_populations; %number of species (populations) active feeders
        C_sp_pass=Cpass_populations;%5; %number of species (population) passive feeders
        param.C_sp_act=C_sp_act;
        param.C_sp_pass=C_sp_pass;
        param.C_sp=C_sp_act+C_sp_pass; %total number of populations
        param.nbr_discr=C_size_classes;
        param.min_cop_act=log10(Cact_min); % minimum size of adult active copepods in log10
        param.max_cop_act=log10(Cact_max); % maximum size of adult active copepods in log10
        param.min_cop_pass=log10(Cpass_min); % minimum size of adult passive copepods in log10
        param.max_cop_pass=log10(Cpass_max); % maximum size of adult passive copepods in log10
        C_grid_act=zeros(param.nbr_discr,C_sp_act); % each column is a population of active copepods, 
                                                    % each row is a size bin, from juveniles (top) to adults (bottom)
        C_grid_pass=zeros(param.nbr_discr,C_sp_pass); %same as above for passive copepods
        Ca_sp_act=logspace(param.min_cop_act,param.max_cop_act,C_sp_act); %vector of adult copepod size
        Ca_sp_pass=logspace(param.min_cop_pass,param.max_cop_pass,C_sp_pass); %vector of adult copepod size
                
        OA_ratio_act=0.01;%0.002; %offspring to adult ratio 0.003 0.0063
        OA_ratio_pass=0.01;%0.02; %offspring to adult ratio 0.01
        for i=1:C_sp_act
            C_grid_act(:,i)=logspace(log10(Ca_sp_act(i).*OA_ratio_act),log10(Ca_sp_act(i)),param.nbr_discr); %copepod grid 0.003
        end
        for i=1:C_sp_pass
            C_grid_pass(:,i)=logspace(log10(Ca_sp_pass(i).*OA_ratio_pass),log10(Ca_sp_pass(i)),param.nbr_discr); %copepod grid 0.003
        end
       
    otherwise
        disp ('wrong size option')
        dbstop
end            

% %cops
C_grid=cat(2,C_grid_act,C_grid_pass);
C_up=C_grid(2:end,:);   %upper sizes
C_dw=C_grid(1:end-1,:); %lower sizes
W=sqrt(C_up.*C_dw);     %Central sizes (geomean);
param.W=cat(1,W,C_up(end,:));       % total grid of copepod weights, last row is the weight of the adult
param.C_up=cat(1,C_up,C_up(end,:)); %we fix it to include the aduts, useful for the feeding kernels afterwards
param.C_dw=cat(1,C_dw,C_up(end,:)); %we fix it to include the aduts, useful for the feeding kernels afterwards
param.Wa=param.W(end,:); %adults mass
param.nbr_cops=length(param.Wa);
param.deltaC=param.C_up-param.C_dw;
param.nbr_stages=length(param.W(:,1));
z=param.C_dw./param.C_up;
param.z=z(:);
param.nbr_Ctot=length(param.W(:));

param.ind_a=param.nbr_stages:param.nbr_stages:param.nbr_stages*param.nbr_cops; %index of adult copepods
param.ind_b=1:param.nbr_stages:param.nbr_stages*param.nbr_cops; %index of newborn nauplii
indexs=1:param.nbr_stages*param.nbr_cops;
param.ind_j=setdiff(indexs,param.ind_a); %indexs of all juveniles stages
param.ind_rest=setdiff(param.ind_j,param.ind_b); %indexs of juveniles stages except newborns
param.ind_act=1:param.nbr_stages*C_sp_act; %indexs of active copepods
param.ind_pass=param.ind_act(end)+1:length(param.W(:)); %index of passive copepods
param.Wvec=param.W(:); % convert the matrix into a column vector



%%%%%%%%%%%%%%%%%%%%%%%%%

%feeding kernels
param.beta_P=500;%100;%500; %Protists
param.sigma_P=1;
param.beta_act=10000;%10000;%1/(10^-3.66); %active feeders
param.sigma_act=1.5;%2;
param.beta_pass=100;%1/10^-2.35; %passive feeders
param.sigma_pass=1;%1.75; 

%this is the function that defines the thershold for mortality by HTL
% no_small=tanh(param.Wvec');
mmax=max(param.W(:));
beta=param.beta_act;        % for plankton and fish at the moment
mshift=max(param.W(:))./beta;
sigma=param.sigma_act;      % for plankton and fish at the moment

%--------------------------------
%passive feeders tau
PL=0.532*param.Wvec.^0.32;
SR=1.801*PL-0.695;
sw=10^0.38.*PL.^0.93;
ratio=SR./sw;
tau=ratio;
tau(ratio<0)=0;
tau(param.ind_act)=1;
param.tau=tau;
%this function is to have a reduced mortality for the passive feeders
p=1/5+tau(param.ind_pass)*(1-1/5); %!!!!!!!!!!!!!!!!!!!
% p(:)=1;

phi=exp(-(log(beta*param.W(:)/mmax)).^2./(2*sigma^2)); %prey  pref function
phi(param.ind_pass)=phi(param.ind_pass)./5; %reduce preference for passives
no_small=phi;%*0.7+0.3; 
no_small(param.Wvec>mshift)=1; %say that any larger than mshift has full preference in mortality htl
param.p=p;
param.no_small=no_small;

param.m_coef=10^-4;%1e-4;5

% % RATES copepods ---------------------------------------------------------- 
I_act=1.37*param.Wvec(param.ind_act).^(-0.25);            %Ingestion rate
k_act=(0.16*param.Wvec(param.ind_act).^(-0.25));          %Metabolic cost
clearance_act=(0.011*param.Wvec(param.ind_act).^(-0.25)); %Clearance rate
d_act=(param.m_coef.*param.Wvec(param.ind_act).^(-0.25)); %mortality rate 0.006

I_pass=0.40*param.Wvec(param.ind_pass).^(-0.25);             %Ingestion rate non-calanoid
k_pass=(0.048*param.Wvec(param.ind_pass).^(-0.25));          %Metabolic cost thomas non-calanoids 
clearance_pass=(0.0052*param.Wvec(param.ind_pass).^(-0.25)); %Clearance rate
d_pass=(param.m_coef.*param.Wvec(param.ind_pass).^(-0.25));  %mortality rate 0.003

tau2=tau;
tau2(param.ind_pass)=0;

param.I=[I_act; I_pass];
param.k=[k_act; k_pass];
param.alpha_c=[clearance_act; clearance_pass];
param.d_c=[d_act; d_pass];


% Efficiencies
param.eff=0.67;
param.rep_eff=0.25;


% PROTIST --------------------------------------------------------------
param.alpha_N=3.75e-5.*param.V.^(-2/3); % nutrient affinity
param.alpha_F=0.0024.*param.V.^(-1/4);  % food affinity
cL=21;                                  % for light affinity
AL=0.000914;                            % for light affinity
param.alpha_L=(AL*param.V.^(2/3).*(1-exp(-cL.*param.V.^(1/3))))./param.V; %light affinity


param.Qcn=5.6; % C:N ratio

V=2.96e7.*param.V.^1.221;
rho=0.024.*V.^1.1;
Q=0.032.*V.^0.76;
mu_inf=4.7.*V.^-0.26;
param.mu_max=mu_inf.*rho./(mu_inf.*Q + rho);
param.mu_max=param.mu_max + param.mu_max.*0.17; % protitst max uptake rate

param.R=0.2.*param.mu_max; %Metabolic cost


param.m=(10^(-3).*param.V.^-0.25)./(param.V_up./param.V_dw); %mortality rate 


% Environment -----------------------------------------------------------
param.kw=0.03;    % attenuation coefficient of water
param.kchl=3e-5;  % attenuation coefficient of Chl plankton
param.remin=0.05; % reminilazitation

% temperature Q10 for physiological rates
param.Q_I=2;   % ingestion rate
param.Q_N=1.5; % nitrogen uptake
param.Q_k=2;   % respiration
param.Tref=15; % reference temperature

%input of individuals "everything is everywhwre"(?)
param.inputP=(0.001*param.V.^-1).*param.delta_V;
param.inputC=(0.001*param.W(1,:).^-1).*param.deltaC(1,:);
param.inputCa=(0.001*0.1*param.W(end,:).^-1).*param.deltaC(end,:);
param.flow=param.mu_max'.*1e-3; % influx rate
param.flowC=param.I(param.ind_b).*1e-3;
param.flowCa=param.I(param.ind_a).*1e-3;

%we asume large copepods to migrate, small no --is not working at the
%momemnt
migration=ones(size(param.W(:)));
migration(param.W>1e1)=1;
param.migration=migration;
param.D=10^-2.5;%0.01;


deltaC=param.deltaC;
deltaC2 = deltaC;
deltaC(end,:)=param.deltaC(param.ind_a-1);
param.deltaC2=deltaC;
deltaratio=param.C_up./param.C_dw;
deltaratio(end,:)=deltaratio(end-1,:);
param.deltaratio=deltaratio;
param.no_small=no_small;

param.Diff_min= 0.2;%  %nutrient Pompei et al: 0.2 
param.remin_frac=0.1;

v=(param.V/(7.6*1e-7)).^(1/0.819);

param.VN=2.9617e+04.*param.V.^(1.1844).*5.6*24*1e-6./param.V;
Pc1=10^(-1.45).*v.^(0.3).*24;
Pc2=10^(-0.44).*v.^(-0.2).*24;
param.VL=min(Pc1,Pc2);

param.VF=10^-0.19.*24./1000.*(param.V./1000).^-0.33;


end