temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')+3;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_plus3            = [temp_winter temp_spring temp_summer temp_autumn temp_annual];


temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')+2;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_plus2            = [temp_winter temp_spring temp_summer temp_autumn temp_annual];

temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')+1;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_plus1            = [temp_winter temp_spring temp_summer temp_autumn temp_annual];

temp1                 = load('meanTemp_slope_MAB_2000_2010.txt');          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_0                = [temp_winter temp_spring temp_summer temp_autumn temp_annual];


temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')-1;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_minus1           = [temp_winter temp_spring temp_summer temp_autumn temp_annual];


temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')-2;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_minus2           = [temp_winter temp_spring temp_summer temp_autumn temp_annual];


temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')-3;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_minus3           = [temp_winter temp_spring temp_summer temp_autumn temp_annual];


temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')-4;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_minus4           = [temp_winter temp_spring temp_summer temp_autumn temp_annual];

temp1                 = load('meanTemp_slope_MAB_2000_2010.txt')-5;          
x                     = temp1-273;
temp_annual           = mean(x);
temp_winter           = (sum(x(337:366))+ sum(x(1:60)))/91;
temp_spring           = mean(x(61:152));
temp_summer           = mean(x(153:244));
temp_autumn           = mean(x(245:336, :));
temp_minus5           = [temp_winter temp_spring temp_summer temp_autumn temp_annual];


%y = {'temp_winter', 'temp_spring', 'temp_summer', 'temp_autumn', 'temp_annual'};
%temp_total            = table(temp_plus3, temp_plus2, temp_plus1, temp_0, temp_minus1, temp_minus2, temp_minus3, temp_minus4, temp_minus5);

temp_total            = [temp_plus3; temp_plus2; temp_plus1; temp_0; temp_minus1; temp_minus2; temp_minus3; temp_minus4; temp_minus5];



% filename = 'temp_SST_2000_2010_slope_mab.xlsx';
% writematrix(temp_total,filename,'Sheet',1,'Range','A1')







