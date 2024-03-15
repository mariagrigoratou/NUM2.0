function [temp,temp_cop, Diff,zmld,Ls,dzdt]=func_seasonal_forcing(param, years)

%load file for the mixed layer depth
% %as in Pompei et al., 2020
%fash=load('mld_fasham_mat.mat');   % 50 lat as in Pompei et al., 2020
% zmld1=fash.alk3; %for Pompei et al., 2020;


% %as in Grigoratou et al., 2023
mld_data              = readtable('slope_mld_2000_2010.xlsx'); %observed mld, Grigoratou et al., 2024 

fash                   = (table2array(mld_data(:,{'mld'}))); % lat- added by MG, May 2021
zmld1=fash;    % for the mld_data Grigoratou


%adapt it to the number of days for the run
zmld1=zmld1(2:end);
zmld1(end)=zmld1(1);
zmld=[zmld1; repmat(zmld1(1:end),years-1,1)]; %mixed layer depth
dzdt=zmld(2:end)-zmld(1:end-1); %change in mixed layer depth
dzdt=[dzdt(1); dzdt]; %correct for vector size

%load temperature file
% %as in Pompei et al., 2020
% temp1=load('tempfile_lat55.txt')+3; % temperature force as in Pompei et al., 2020
% temp1=temp1(2:end);
% temp=[temp1(1); repmat(temp1(1:end),years,1)];
% temp=temp-273;

% %as in Grigoratou et al., 2023
temp1=load('meanTemp_slope_MAB_2000_2010.txt');% 365 days
temp=[temp1(1); repmat(temp1(1:end),years,1)]; 

% Heatwave parameters 
temp_cop = temp; %heatwave for either only protists or copepods
warming = 4.0; %Heatwave temp


% %HEATWAVES both cop and protists -- is being called in the z_ode_copepod_model for the T_inn_kelvin
% temp(18219:18250) = temp(18219:18250) + warming; % December 50th year  
% temp(18251:18281) = temp(18251:18281) + warming; % January 
% temp(18282:18309) = temp(18282:18309) + warming; % February
% temp(18310:18340) = temp(18310:18340) + warming; % March
% temp(18341:18370) = temp(18341:18370) + warming; % April
% temp(18371:18401) = temp(18371:18401) + warming; % May
% temp(18402:18431) = temp(18402:18431) + warming; % June
% temp(18432:18462) = temp(18432:18462) + warming; % July
% temp(18463:18493) = temp(18463:18493) + warming; % August
% temp(18494:18523) = temp(18494:18523) + warming; % September
% temp(18524:18554) = temp(18524:18554) + warming; % October
% temp(18555:18584) = temp(18555:18584) + warming; % November


% %HW cops %% check the HWs only -- is being called in the z_ode_copepod_model for the T_inn_kelvin_HW 
%  temp_cop(18219:18250) = temp_cop(18219:18250) + warming; % December 50th year
%  temp_cop(18251:18281) = temp_cop(18251:18281) + warming; % January 
%  temp_cop(18282:18309) = temp_cop(18282:18309) + warming; % February
%  temp_cop(18310:18340) = temp_cop(18310:18340) + warming; % March
%  temp_cop(18341:18370) = temp_cop(18341:18370) + warming; % April
%  temp_cop(18371:18401) = temp_cop(18371:18401) + warming; % May
%  temp_cop(18402:18431) = temp_cop(18402:18431) + warming; % June
%  temp_cop(18432:18462) = temp_cop(18432:18462) + warming; % July
%  temp_cop(18463:18493) = temp_cop(18463:18493) + warming; % August
%  temp_cop(18494:18523) = temp_cop(18494:18523) + warming; % September
%  temp_cop(18524:18554) = temp_cop(18524:18554) + warming; % October
%  temp_cop(18555:18584) = temp_cop(18555:18584) + warming; % November


%mixing rate (input of N)
Diff=(param.Diff_min+max(0,dzdt))./zmld;

%params light function
clouds=1-(5/8);%cloud cover % Pompei default: 1-(6/8)

kyear=0;
day_type=1;
tsol=1:1:365;
Lss=daily_insolation(kyear,40,tsol,day_type); %surface light
Ls=Lss.*clouds.*0.4.*4.6; % as in Pompei et al 2020 . correct for cloud cover, convert to PAR (PAR/I=0.4) and convert to microeinsteins (1Wm^-2=4.6microeinteins)


%correct format - comment by Pompei 
Ls=Ls';
%Ls=[Ls(1); repmat(Ls(1:end),years,1)]; %as in Pompei et al 2020
Ls= repmat(Ls(1:end),years,1); %remove the 366th day - added by MG Sep 2021
Ls=Ls';

%added my MG, May 2021
% figure
% plot(Lss)
% hold on 
% plot(Ls)
% % hold on 
% % plot(Ls1)
% % hold on 
% % plot(Ls2)
% legend ('no clouds', 'clouds 1 -5/8');
% title (' Light (in microeinsteins)in 40 Lat, with different clouds ');
% hold off


end


