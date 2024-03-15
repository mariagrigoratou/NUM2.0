% mixed layer as in Pompei et al 2020
function dydt=mixed_layer(t,y,H1)
h=y(1);
T=y(2);

if t<1
H1=0;
else
H1=H1(ceil(t));
end
% H1=daily_insolation(param.kyear,param.lat,t,param.day_type)*(1/8);
% H1= light_brock_ocean(t,param.lat,1);
% H1=light_brock_ocean(t,param.lat,1);
H=H1* 86400;


rho=1035;
rhoa=1.2466;
wind=10;
alpha=200e-6;
g=9.82;
C=4200;
me=1.5e-3;
T0=4+273;
deltaT = T - T0;

cd = (0.8+0.0014*wind)*1e-3; % Drag coefficient
m=me*(rho/(rhoa*cd))^0.5;
u = sqrt(cd*rhoa*wind^2/rho);%friction velocity

% turb=((86400 *2*m*u^3)/(g*alpha)); %turbulent kinetic energy

test=(((86400 *2*m*u^3)/(g*alpha))-((H*h)/(rho*C)))/(deltaT*h);

if test>0
    dhdt=test;
    dTdt=(2/h^2)*((H*h/(rho*C))-((86400*2*m*u^3)/(g*alpha)));
    
else   
    h_t=(2 * 86400*m*u^3)/((g*alpha*H)/(rho*C));
    dhdt = h_t - h;
    
    dTdt=0.5*((H/(rho*C))^2)/((86400 *2*m*u^3)/(g*alpha));
    
end

dydt=[dhdt; dTdt];
end