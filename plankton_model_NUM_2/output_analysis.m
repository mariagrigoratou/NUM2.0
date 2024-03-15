%
C_lastyr = C(17885:18250, :);
P_lastyr = P(17885:18250, :);
Cadult_lastyr = C(17885:18250, 8:8:704);
writematrix(C_lastyr, 'C_slope_mab_mld_1999_average_sst_2000_2010_minus5degrees.csv')
writematrix(P_lastyr, 'P_slope_mab_mld_1999_average_sst_2000_2010_minus5degrees.csv')
writematrix(Cadult_lastyr, 'P_slope_mab_mld_1999_average_sst_2000_2010_minus5degrees.csv')
% %% call variables from multiple .mat 
% 
%md1                 = load('model_output\model_output_1');
%md2                 = load('model_output\model_output_2');
% md3                 = load('model_output\model_output_3');
% md4                 = load('model_output\model_output_4');
% md5                 = load('model_output\model_output_5');
% md6                 = load('model_output\model_output_6');
% md7                 = load('model_output\model_output_7');
% md8                 = load('model_output\model_output_8');
% md9                 = load('model_output\model_output_9');
% md10                = load('model_output\model_output_10');
% md11                = load('model_output\model_output_11');
% md12                = load('model_output\model_output_12');
% md13                = load('model_output\model_output_13');
% md14                = load('model_output\model_output_14');
% md15                = load('model_output\model_output_15');
% md16                = load('model_output\model_output_16');
% md17                = load('model_output\model_output_17');
% md18                = load('model_output\model_output_18');
% md19                = load('model_output\model_output_19');
% md20                = load('model_output\model_output_20');
%% plots
%C = cat(1, md1.C, md2.C);%, md3.C, md4.C, md5.C, md6.C, md7.C, md8.C, md9.C, md10.C);
%C2 = cat(1, md11.C, md12.C, md13.C, md14.C,md15.C, md16.C, md17.C, md18.C, md19.C, md20.C);
%P = cat(1, md1.P, md2.P);%, md2.P, md3.P, md4.P, md5.P, md6.P, md7.P, md8.P, md9.P, md10.P);
%P2 = cat(1, md11.P, md12.P, md13.P, md14.P,md15.P, md16.P, md17.P, md18.P, md19.P, md20.P);
%P = cat(1, P1, P2);
%C = cat(1, C1,C2);

%%
time              = 17885:18250;                      %in days (first 5 years: 1:1825, 5 last years:16425:18250 time, 50 yrs: 1:18250
% actSpecies        = (1:8:519);               % 1-512:active cops, 513-704: passive cops 
% paSpecies         = (520:8:704);
% stages            = length(species(:));
% warm_species      = 184;                       % warm species
% dydt_species      = 113 + species;             % 1: nutrient, 2-113:prot, 114-625: active cops, 626- 817: passive cop 
% dydt_warm_species = 113 + warm_species;        % 1: nutrient, 2-113:prot, 114-625: active cops, 626- 817: passive cop 
% gamma_species     = species-1;                 % species for juv, species-1 for adult
% repro             = round(stages/8.0);        % repro rates, 1-64: active cops, 65- 88: passive cops works for seeding as well
% gamma_warm_species= warm_species-1;            % species for juv, species-1 for adult
% warm_repro        = round(warm_species/8.0);   % repro rates, 1-64: active cops, 65- 88: passive cops works for seeding as well
aSpec0            = (8:8:71);
aSpec4            = (72:8:135);
aSpec8            = (136:8:199);
aSpec12           = (200:8:263);
aSpec16           = (264:8:327);
aSpec20           = (328:8:391);
aSpec24           = (392:8:455);
aSpec28           = (456:8:519);

pSpec0            = (520:8:543);
pSpec4            = (544:8:567);
pSpec8            = (568:8:591);
pSpec12           = (592:8:615);
pSpec16           = (616:8:639);
pSpec20           = (640:8:663);
pSpec24           = (664:8:687);
pSpec28           = (688:8:704);

pro0               = (1:15);
pro4              = (16:29);
pro8              = (30:43);
pro12             = (44:57);
pro16             = (58:71);
pro20             = (72:85);
pro24             = (86:99);
pro28             = (100:112);
% figure
% %subplot(2,4,1)
% plot(C(time,species))
% hold on
% plot(C8.C(time,warm_species))
% legend ('8 temp', '8oC')
% title('C');
% hold off
% data_active              = readtable('C:\Users\mgrigoratou\Documents\GitHub\Copepod_sizebased_model-master_no_fp\plankton_model_v2\mld_1997\\mld_shlef_passive_1999.xlsx');
% 
% figure  
% scatter(data_active.doy,log10(data_active.max_bio_mld*0.001))
% hold on 
% scatter(data_active.doy,log10(data_active.min_bio_mld*0.001) ,'*')
% legend ('max', 'min')
% title ('log10biomass Oithona (mgC m3) 1999')
% hold off
%% active cops

figure
subplot(2,4,1)
plot(C(time,aSpec0))
%area(C(time,aSpec0))
legend
title('active cops 0oC');
xlim([0 366]);
% legf=legend(plots([1 2 3 4 5 6]),flip({'m < 10^{-2} ','10^{-2} \leq m < 10^{-1}','10^{-1} \leq m < 10^{0}','10^{0} \leq m < 10^{1}','10^{1} \leq m < 10^{2}','10^{2} \leq m'}));
% legend boxoff
% legf.FontSize = 8;
hold off

subplot(2,4,2)
plot(C(time,aSpec4))
legend
title('active cops 4oC');
xlim([0 366]);
hold off

subplot(2,4,3)
plot(C(time,aSpec8))
% hold on
% plot(cpr3.doy, cpr3.min_biomass, '*')
legend
title('active cops 8oC');
xlim([0 366]);
hold off

subplot(2,4,4)
plot(C(time,aSpec12))
legend
title('active cops 12oC');
% hold on
% plot(cpr3.doy, cpr3.min_biomass, '*')
xlim([0 366]);
hold off

subplot(2,4,5)
plot(C(time,aSpec16))
legend
title('active cops 16oC');
xlim([0 366]);
hold off

subplot(2,4,6)
plot(C(time,(aSpec20)))
%plot(C(time,(640:663)))
legend
% hold on 
% scatter(data_active.doy,(data_active.min_bio_mld*0.001))
% hold on 
% scatter(data_active.doy,(data_active.max_bio_mld*0.001))
title('active cops 20oC ');
xlim([0 366]);
hold off

subplot(2,4,7)
plot(C(time,aSpec24))
legend
title('active cops 24oC');
xlim([0 366]);
hold off

subplot(2,4,8)
plot(C(time,aSpec28))
legend
title('active cops 28oC');
xlim([0 366]);
hold off

%% passive cops

figure
subplot(2,4,1)
plot(C(time,pSpec0))
legend
title('passive cops 0oC');
xlim([0 366]);
hold off

subplot(2,4,2)
plot(C(time,pSpec4))
legend
title('passive cops 4oC');
xlim([0 366]);
hold off

subplot(2,4,3)
plot(C(time,pSpec8))
legend
title('passive cops 8oC');
xlim([0 366]);
hold off

subplot(2,4,4)
plot(C(time,pSpec12))
legend
title('passive cops 12oC');
xlim([0 366]);
hold off

subplot(2,4,5)
plot(C(time,pSpec16))
legend
title('passive cops 16oC');
xlim([0 366]);
hold off

subplot(2,4,6)
plot(C(time,pSpec20))
legend
title('passive cops 20oC');
xlim([0 366]);
hold off

subplot(2,4,7)
plot(C(time,pSpec24))
legend
title('passive cops 24oC');
xlim([0 366]);
hold off

subplot(2,4,8)
plot(C(time,pSpec28))
legend
title('passive cops 28oC');
xlim([0 366]);
hold off


%% protists and nutrients
time              = 17885:18250;     
figure
subplot(2,4,1)
plot(P(time,pro0))
%legend
title('Protists 0oC');
xlim([0 366]);
hold off

subplot(2,4,2)
plot(P(time,pro4))
%legend
title('Protists 4oC');
xlim([0 366]);
hold off

subplot(2,4,3)
plot(P(time,pro8))
%legend
title('Protists 8oC');
xlim([0 366]);
hold off

subplot(2,4,4)
plot(P(time,pro12))
%legend
title('Protists 12oC');
xlim([0 366]);
hold off

subplot(2,4,5)
plot(P(time,pro16))
%legend
title('Protists 16oC');
xlim([0 366]);
hold off

subplot(2,4,6)
plot(P(time,pro20))
%legend
title('Protists 20oC');
xlim([0 366]);
hold off

subplot(2,4,7)
plot(P(time,pro24))
%legend
title('Protists 24oC');
xlim([0 366]);
hold off

subplot(2,4,8)
plot(P(time,pro28))
%legend
title('Protists 28oC');
xlim([0 366]);
hold off


figure
plot(P(time,:))
title ('Protists')'
xlim([0 366]);

figure
plot(N(time,:))
title ('Nutrients')'
xlim([0 366]);

figure
plot(T(time,:))
title ('Temp')'
xlim([0 366]);

% subplot(1,2,1)
% plot(C(time,aSpec8))
% hold on
% plot(C(time,pSpec8))
% hold on
% plot(cpr3.doy, cpr3.min_biomass, '*')
% legend
% title('active cops 8oC');
% ylim([1e-06 6e-04])
% hold off
% 
% subplot(1,2,2)
% plot(C(time,pSpec12))
% legend
% title('passive cops 12oC');
% hold on
% plot(C(time,pSpec12))
% hold on
% plot(cpr3.doy, cpr3.min_biomass, '*')
% ylim([1e-06 6e-04])
% hold off


