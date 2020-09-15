%John Furumo
%ENAE 601 section 0101
%Final Project
%due 12/2/19

clear;
clf;
clc;
close all;

format longg;
%format short;
radius_earth=6378; %km
G=6.67428E-020; %km^3/kg*s^2
Spacecraft_mass=100; %4x 12U cubesats of 25kg each
n=10; %number of bodies: sun + 9 planets


%solar system body
body{1}=('Sun.txt');
body{2}=('Mercury.txt');
body{3}=('Venus.txt');
body{4}=('Earth.txt');
body{5}=('Mars.txt');
body{6}=('Jupiter.txt');
body{7}=('Saturn.txt');
body{8}=('Uranus.txt');
body{9}=('Neptune.txt');
body{10}=('Pluto.txt');

%solar system mu values (km^3/s^2)
mu(1)=132712440041.93938; %Sun
mu(2)=22031.78;           %Mercury
mu(3)=324858.592;         %Venus
mu(4)=398600.435436;      %Earth
mu(5)=42828.375214;       %Mars
mu(6)=126686534.911;      %Jupiter
mu(7)=37931207.8;         %Saturn
mu(8)=5793951.322;        %Uranus
mu(9)=6835099.5;          %Neptune
mu(10)=869.33907803;      %Pluto
%switch to column vector
mu=mu';

for p=1:n
%folder pathway
folder=('C:\Users\jfuru\Desktop\UMD\Fall 2019\ENAE 601\final project\Solar System Ephemeris\2020-2034\');
%full pathway
path=strcat(folder,body{p});
%open TLE file
fID = fopen(path,'r');
%read data from ephemeris file line by line
%find end of header line
%find ephermeris data start line
k=1;
while feof(fID) == 0
L{k} = fgetl(fID);
     if strfind(L{k},'$$SOE') ~= 0
        data_start=k;
     end
     if strfind(L{k},'$$EOE') ~= 0
        data_end=k;
     end 
k=k+1;
end
%read state vector
comp_length=22;
for j=1:2
    for i=1:3
        R(i+3*(j-1),:)=str2num(L{(data_start+1+j)}(((5*i)+(comp_length-1)*(i-1)):(5+(comp_length-1))*i));
    end
end
SS{p,:}={body{p}(1:end-4),R'};
fclose(fID);
end

%number of ephemeris data points
n_ephemeris=str2num(L{(data_end-3)}(1:7))-str2num(L{(data_start+1)}(1:7));
%celldisp(SS)

%Cubesat mission to Saturn/Titan
%start date: January 1, 2020
%end date:   January 1, 2034 (5114 days later)
G=6.67428E-020; %km^3/kg*s^2
tf=n_ephemeris*86400; %s
%one state vector per day
tspan=[0:86400:tf]; %s
%tspan=[0:1:n_ephemeris]; %day
tolerance=1e-013;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
n=10; %number of bodies: sun + 9 planets

%Solar System initial state vector
for i=1:n
SS_state(1,(6*(i-1)+1):(6*i))=cell2mat((SS{i}(2)));
end
SS_state=SS_state';

%Numerical integration of solar system
[t,system_N]=ode45(@integrate_NBP,tspan,SS_state,options,n,mu);

%3D plot of solar system
figure('Name','Solar System');
plot3(system_N(:,1),system_N(:,2),system_N(:,3),'y*','linewidth',2) %sun
hold on;
plot3(system_N(:,7),system_N(:,8),system_N(:,9),'r','linewidth',2) %mercury
plot3(system_N(:,13),system_N(:,14),system_N(:,15),'y','linewidth',2) %venus
plot3(system_N(:,19),system_N(:,20),system_N(:,21),'b','linewidth',2) %earth
plot3(system_N(:,25),system_N(:,26),system_N(:,27),'r','linewidth',2) %mars
plot3(system_N(:,31),system_N(:,32),system_N(:,33),'color',[0.8500, 0.3250, 0.0980],'linewidth',2) %jupiter
plot3(system_N(:,37),system_N(:,38),system_N(:,39),'color',[0.9290, 0.6940, 0.1250],'linewidth',2) %saturn
%plot3(system_N(:,43),system_N(:,44),system_N(:,45),'g','linewidth',2) %uranus
%plot3(system_N(:,49),system_N(:,50),system_N(:,51),'b','linewidth',2) %neptune
%plot3(system_N(:,55),system_N(:,56),system_N(:,57),'k','linewidth',2) %pluto
%plot Earth initial and final positions
plot3(system_N(1,19),system_N(1,20),system_N(1,21),'go','linewidth',2,'markersize',10)
plot3(system_N(end,19),system_N(end,20),system_N(end,21),'rx','linewidth',2,'markersize',10)
%plot Saturn initial and final positions
plot3(system_N(1,37),system_N(1,38),system_N(1,37),'go','linewidth',2,'markersize',10)
plot3(system_N(end,37),system_N(end,38),system_N(end,37),'rx','linewidth',2,'markersize',10)
axis equal;
title({'Solar System - Heliocentric Inertial Coordinate Frame','NBP Simulation 01/01/2020 - 01/01/2034'})
legend([{'Sun'},{'Mercury'},{'Venus'},{'Earth'},{'Mars'},{'Jupiter'},{'Saturn'},{'Earth initial position'},{'Earth final position'},{'Saturn initial position'},{'Saturn final position'}])
xlabel('X distance (km)')
ylabel('Y distance (km)')
zlabel('Z distance (km)')
xlim([-1.5e+009 1.5e+009]);
ylim([-1.5e+009 1.5e+009]);
az = 0;
el = 90;
view(az, el);

%3D plot of Earth, Saturn/Titan system
figure('Name','Solar System');
%plot Sun
p1=plot3(system_N(:,1),system_N(:,2),system_N(:,3),'yo','linewidth',2);
hold on;
%plot Earth
p2=plot3(SS_state(19),SS_state(20),SS_state(21),'go','linewidth',2);
p3=plot3(system_N(:,19),system_N(:,20),system_N(:,21),'b');
p4=plot3(system_N(end,19),system_N(end,20),system_N(end,21),'ro','linewidth',2);
%plot Saturn
p5=plot3(SS_state(37),SS_state(38),SS_state(39),'go','linewidth',2);
p6=plot3(system_N(:,37),system_N(:,38),system_N(:,39),'y');
p7=plot3(system_N(end,37),system_N(end,38),system_N(end,39),'ro','linewidth',2);
%plot3(system_N(:,61),system_N(:,62),system_N(:,63),'k')
axis equal;
title({'Solar System','NBP Simulation 01/01/2020 - 01/01/2026'})
legend([p1 p2 p4 p5 p7],{'Sun','Earth Initial Location','Earth Final Location','Saturn Initial Location','Saturn Final Location'})
xlabel('X distance (km)')
ylabel('Y distance (km)')
zlabel('Z distance (km)')
% az = 0;
% el = 90;
% view(az, el);

%distance from Earth to Saturn
for i=1:length(t)
    distance(i)=sqrt((system_N(i,19)-system_N(i,37))^2+(system_N(i,20)-system_N(i,38))^2+(system_N(i,21)-system_N(i,39))^2);
end
distance=distance';
%2D plot of Earth - Saturn/Titan distance
figure('Name','Earth to Saturn/Titan');
plot(t,distance)
title({'Distance from Earth to Saturn/Titan'})
xlabel('Time (seconds)')
ylabel('Distance (km)')
%xlim([-1.5e+009 1.5e+009]);
ylim([0 2e+009]);
%%
%Patched Conics Method
%mu_sun = mu(1)
%mu_earth = mu(4)
R_earth = 149.6e+06; %mean orbital radius (km)
perihelion_E = 147.1e+06; %km
aphelion_E = 152.1e+06; %km
a_E=(perihelion_E+aphelion_E)/2; %km
e_E=aphelion_E/a_E-1
nE = 0.9856474; %mean daily motion (deg/day)
nE = nE*(2*pi)/(360*86400); %mean motion (rad/s)
orbitE=R_earth*ones(1000,1); %km
%mu_saturn = mu(7)
R_saturn = 1433.55e+06; %mean orbital radius (km)
perihelion_S = 1352.6e+06; %km
aphelion_S = 1514.5e+06; %km
a_S=(perihelion_S+aphelion_S)/2; %km
e_S=aphelion_S/a_S-1
nS = 0.0334979; %mean daily motion (deg/day)
nS = nS*(2*pi)/(360*86400); %mean motion (rad/s)
orbitS=R_saturn*ones(1000,1); %km
%Earth - Saturn synodic period
Ts=(2*pi)/abs(nE-nS); %s
Ts_day=Ts/86400; %days
%Heliocentric Hohmann Transfer from Earth to Saturn
aT=(R_earth+R_saturn)/2; %km
eT=1-R_earth/aT;
ToF=(pi*sqrt(aT^3/mu(1))); %s
r1=system_N(1,19:21);
r2=system_N(1,37:39);
delta_theta_0=acos(dot(r1,r2)/(norm(r1)*norm(r2))); %rad
zeta=nS*ToF; %rad/s
delta_theta=zeta-pi; %rad
%wait time
t_wait=(delta_theta-delta_theta_0+2*pi)/abs(nE-nS);
%transit windows
win_d=14*365-ToF/(86400); %window (days)
n_opp=floor((win_d-t_wait/86400)/Ts_day);
%first transit opportunity
nu0=0;
[t_wait,dv1,dv2,dv_rendezvous,ToF,T_syn] = Rendezvous_Maneuver(R_earth,R_saturn,delta_theta_0,nu0,1,mu(1))
nu0=nE*(t_wait+Ts);
%subsequent transit opportunities
% for opp=2:n_transit_opp
% 
% 
% [t_wait,~,~,~] = Rendezvous_Maneuver(R_earth,R_saturn,delta_theta,nu0,opp,mu(1))
% nu0=nE*Ts;
% 
% end

%%

% %Earth orbital velocity (circular orbit assumption)
vE=sqrt(mu(1)/R_earth) %km/s
%Spacecraft transfer ellipse velocity at periapsis (Earth)
vTE=sqrt(((2*mu(1))/R_earth)-(2*mu(1))/(R_earth+R_saturn)) %km/s
%difference between spacecraft and earth velocity
v_inf_E=vTE-vE %km/s
%Saturn orbital velocity (circular orbit assumption)
vS=sqrt(mu(1)/R_saturn) %km/s
%Spacecraft transfer ellipse velocity at apoapsis (Earth)
vTS=sqrt(((2*mu(1))/R_saturn)-(2*mu(1))/(R_earth+R_saturn)) %km/s
%difference between Saturn and spacecraft velocity
v_inf_S=vS-vTS %km/s

%circular Earth parking orbit
r1=7000; %km
v1=sqrt(mu(4)/r1) %km/s
%hyperbolic Earth escape orbit
vH1=sqrt(v_inf_E^2+(2*mu(4))/r1) %km/s
%Earth departure delta V
dV1=vH1-v1 %km/s
%Earth departure delta-v (601 formula)
C3=v_inf_E^2
dv_dep=sqrt(C3+2*mu(4)/r1)-sqrt(mu(4)/r1)
%Earth sphere of influence
r_SOI_earth=R_earth*(mu(4)/mu(1))^(2/5) %km
%plot Earth SOI departure
e_dep=1+(r1*v_inf_E^2)/mu(4)
beta_dep=acosd(1/e_dep)
% figure('name','Earth SOI Departure')
% %plot3(0,0,0,'bo')
% nu=linspace(0,2*pi,1000);
% nu=nu'; %rad
% polarplot(0,0,'LineWidth',2,'color','g')
% hold on
% polarplot(nu,radius_earth*ones(1000,1),'b-','LineWidth',2)
% polarplot(nu,r1*ones(1000,1),'LineWidth',2,'color','r')
% polarplot(nu,r_SOI_earth*ones(1000,1),'c--','LineWidth',2)



%circular Saturn parking orbit
r2=70000; %km
v2=sqrt(mu(7)/r2) %km/s
%hyperbolic Saturn capture orbit
vH2=sqrt(v_inf_S^2+(2*mu(7))/r2) %km/s
%Saturn capture delta V
dV2=vH2-v2 %km/s
%Saturn arrival delta-v (601 formula)
C3=v_inf_S^2
dv_arr=sqrt(C3+2*mu(7)/r2)-sqrt(mu(7)/r2)
%Saturn sphere of influence
r_SOI_saturn=R_saturn*(mu(7)/mu(1))^(2/5) %km
%total delta V
dV_tot=dV1+dV2 %km/s

%propellant mass required (work backwards through maneuvers)
Isp=300; %s
g0=9.81; %m/s^2
%Saturn arrival maneuvers
m_i=Spacecraft_mass*exp((1000*dv_arr)/(Isp*g0)); %kg

%Solar System Planet Circular Orbit Approximation
% figure('name','Planetary Orbits')
% nu=linspace(0,2*pi,1000);
% nu=nu'; %rad
% Rad_earth=R_earth*ones(1000,1); %km
% Rad_saturn=R_saturn*ones(1000,1); %km
% %circular orbit approximations
% polarplot(nu,Rad_earth,'LineWidth',2,'color','g','DisplayName','Earth Orbit (circular)');
% hold on
% polarplot(nu,Rad_saturn,'LineWidth',2,'color','g','DisplayName','Saturn Orbit (circular)');
% %elliptical orbits
% aE=(perihelion_E+aphelion_E)/2; %km
% eE=1-perihelion_E/aE; %unitless
% aS=(perihelion_S+aphelion_S)/2; %km
% eS=1-perihelion_S/aS; %unitless
% for i=1:1000
%     Earth_orbit_elliptical(i)=aE*(1-eE^2)/(1+eE*cos(nu(i)));
%     Saturn_orbit_elliptical(i)=aS*(1-eS^2)/(1+eS*cos(nu(i)));
% end
% polarplot(nu,Earth_orbit_elliptical,'LineWidth',2,'color','b','DisplayName','Earth Orbit (elliptical)');
% polarplot(nu,Saturn_orbit_elliptical,'LineWidth',2,'color','b','DisplayName','Saturn Orbit (elliptical)');
% legend()


%% Cubesat mission to Saturn/Titan
%start date: January 1, 2020
%end date:   January 1, 2026 (2192 days later)
G=6.67428E-020; %km^3/kg*s^2
m_SC=100; %4x 12U cubesats of 25kg each
%radius of Earth parking orbit
r_park=42000; %km


%% Hohmann Transfer Opportunities
%Plot Solar System Initial State
% figure('name','Solar System')
% %subplot(4,2,1);
% form=[{'*'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'}];
% %color=[[1 1 0],[0.6350 0.0780 0.1840],{[0.8500 0.3250 0.0980]},{[0 0.4470 0.7410]},{[1 0 0]},{[0.8500, 0.3250, 0.0980]},{[0.9290, 0.6940, 0.1250]},{[0.4660 0.6740 0.1880]},{[0.3010 0.7450 0.9330]},{[0 0 0]},{[0 0 0]}];
% body=[{'Sun'},{'Mercury'},{'Venus'},{'Earth'},{'Mars'},{'Jupiter'},{'Saturn'},{'Uranus'},{'Neptune'},{'Pluto'},{'spacecraft'}];
% NBP_plot(n,system_N,'Solar System Body Motion 2020 - 2034',form)



for opp=1:n_opp
%one state vector per day
t0=round(t_wait+(opp-1)*T_syn); %s
tf=round(t_wait+(opp-1)*T_syn+ToF); %s
tspan=[t0:86400:tf]; %s
%tspan=[0:1:n_ephemeris]; %day
tolerance=1e-006;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
n=11; %number of bodies: sun + 9 planets + spacecraft
mu(11)=m_SC*G;
%Solar System initial state vector
[val,t_row(opp)]=min(abs(t-t0));
r_E(opp,:)=system_N(t_row(opp),19:21);
v_hat(opp,:)=system_N(t_row(opp),22:24)./norm(system_N(t_row(opp),22:24));  
r_hat(opp,:)=system_N(t_row(opp),19:21)./norm(system_N(t_row(opp),19:21));
SS_state_SC=system_N(t_row(opp),:);
SS_state_SC(61:63)=SS_state_SC(19:21)+((r_SOI_earth)*r_hat(opp,:));
SS_state_SC(64:66)=SS_state_SC(22:24)+((dv1)*v_hat(opp,:));

SS_state_PC=SS_state_SC(61:66);
%Numerical integration of solar system
[t_SC,system_N_SC]=ode45(@integrate_NBP,tspan,SS_state_SC,options,n,mu);

%Patched Conics 2BP Integration
[t_PC,system_N_PC]=ode45(@integrate_2BP,tspan,SS_state_PC,options,mu(1));

%3D plot of NBP system
figure('Name','Solar System');
plot3(system_N(:,1),system_N(:,2),system_N(:,3),'y*','linewidth',2) %sun
hold on;
plot3(system_N(:,7),system_N(:,8),system_N(:,9),'r','linewidth',2) %mercury
plot3(system_N(:,13),system_N(:,14),system_N(:,15),'y','linewidth',2) %venus
plot3(system_N(:,19),system_N(:,20),system_N(:,21),'b','linewidth',2) %earth
plot3(system_N(:,25),system_N(:,26),system_N(:,27),'r','linewidth',2) %mars
plot3(system_N(:,31),system_N(:,32),system_N(:,33),'color',[0.8500, 0.3250, 0.0980],'linewidth',2) %jupiter
plot3(system_N(:,37),system_N(:,38),system_N(:,39),'color',[0.9290, 0.6940, 0.1250],'linewidth',2) %saturn
axis equal;
% title({'Solar System - Heliocentric Inertial Coordinate Frame','NBP Simulation 01/01/2020 - 01/01/2034'})
% legend([{'Sun'},{'Mercury'},{'Venus'},{'Earth'},{'Mars'},{'Jupiter'},{'Saturn'},{'Earth initial position'},{'Earth final position'},{'Saturn initial position'},{'Saturn final position'}])
name=char(strcat('Transfer Opportunity',{' '},num2str(opp)));
title(name)
xlabel('X distance (km)')
ylabel('Y distance (km)')
zlabel('Z distance (km)')
xlim([-1.5e+009 1.5e+009]);
ylim([-1.5e+009 1.5e+009]);
az = 0;
el = 90;
view(az, el);
%figure('name','Rendezvous Opportunity')
%subplot(4,2,opp+1);
%form=[{'o'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'}];
%color=[{'y'},{'r'},{'y'},{'b'},{'r'},{'y'},{'y'},{'g'},{'c'},{'k'},{'k'}];
%body=[{'Sun'},{'Mercury'},{'Venus'},{'Earth'},{'Mars'},{'Jupiter'},{'Saturn'},{'Uranus'},{'Neptune'},{'Pluto'},{'spacecraft'}];
plot3(system_N_SC(:,61),system_N_SC(:,62),system_N_SC(:,63),'c-','linewidth',2)
plot3(system_N_PC(:,1),system_N_PC(:,2),system_N_PC(:,3),'g-','linewidth',2)
plot3(system_N_SC(1,19),system_N_SC(1,20),system_N_SC(1,21),'o','linewidth',2,'markersize',10)
plot3(system_N_SC(end,37),system_N_SC(end,38),system_N_SC(end,39),'x','linewidth',2,'markersize',10)
theta=linspace(0,2*pi,1000)
plot3(R_saturn*cos(theta),R_saturn*sin(theta),zeros(1,1000),'y')
%HT_error(opp,1)=norm(system_N_SC(end,37:39)-system_N_SC(end,61:63));
%HT_error(opp,2)=norm(system_N_SC(end,37:39)-system_N_PC(end,1:3));
end
%circular appoximation of Saturn orbit
theta=linspace(0,2*pi,1000)
plot3(R_saturn*cos(theta),R_saturn*sin(theta),zeros(1,1000),'y')
%tangential velocity vector
[~,p_v_hat1] = plot_vector(100000000*v_hat1,r_E1,'r','--','x')
% % 
% legend([{'Sun'},{'Mercury'},{'Venus'},{'Earth'},{'Mars'},{'Jupiter'},{'Saturn'},{'Earth Departure'},{'Saturn Arrival'},{'Lambert Trajectory (long way)'},{'Lambert Trajectory (short way)'}]);

%% Lambert solver
%ToF=192057892.980357; %s (Hohmann Transfer)
%tof=2920; %days
%ToF=tof*86400;
tof=round(ToF/86400); %days
day0=round(t_wait/86400);
%day0=2190;
r1=system_N(day0,19:21);
r2=system_N(day0+tof,37:39);
[v1l,v2l] = Lambert_solver(r1,r2,ToF,'long',mu(1));
state_Lam_l=[r1';v1l'];
Lam_l=cart2oe(state_Lam_l,mu(1));
[v1s,v2s] = Lambert_solver(r1,r2,ToF,'short',mu(1));
state_Lam_s=[r1';v1s'];
Lam_s=cart2oe(state_Lam_s,mu(1));

%3D plot of NBP system
figure('Name','Solar System');
plot3(system_N(:,1),system_N(:,2),system_N(:,3),'y*','linewidth',2) %sun
hold on;
plot3(system_N(:,7),system_N(:,8),system_N(:,9),'r','linewidth',2) %mercury
plot3(system_N(:,13),system_N(:,14),system_N(:,15),'y','linewidth',2) %venus
plot3(system_N(:,19),system_N(:,20),system_N(:,21),'b','linewidth',2) %earth
plot3(system_N(:,25),system_N(:,26),system_N(:,27),'r','linewidth',2) %mars
plot3(system_N(:,31),system_N(:,32),system_N(:,33),'color',[0.8500, 0.3250, 0.0980],'linewidth',2) %jupiter
plot3(system_N(:,37),system_N(:,38),system_N(:,39),'color',[0.9290, 0.6940, 0.1250],'linewidth',2) %saturn
plot3(system_N(day0,19),system_N(day0,20),system_N(day0,21),'o','linewidth',2,'markersize',10)
plot3(system_N(day0+tof,37),system_N(day0+tof,38),system_N(day0+tof,39),'x','linewidth',2,'markersize',10)
%Lambert propagation
tspan=[0:86400:ToF]; %day
tolerance=1e-006;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
[t_Lam_l,system_N_Lam_l]=ode45(@integrate_2BP,tspan,state_Lam_l,options,mu(1));
[t_Lam_s,system_N_Lam_s]=ode45(@integrate_2BP,tspan,state_Lam_s,options,mu(1));
plot3(system_N_Lam_l(:,1),system_N_Lam_l(:,2),system_N_Lam_l(:,3),'c-','linewidth',2)
plot3(system_N_Lam_s(:,1),system_N_Lam_s(:,2),system_N_Lam_s(:,3),'g-','linewidth',2)
theta=linspace(0,2*pi,1000)
plot3(R_saturn*cos(theta),R_saturn*sin(theta),zeros(1,1000),'y')
%Lambert delta-v
dv_l_1=norm(system_N(day0,22:24)-v1l);
dv_s_1=norm(system_N(day0,22:24)-v1s);
dv_l_2=norm(system_N(day0+tof,40:42)-v2l);
dv_s_2=norm(system_N(day0+tof,40:42)-v2s);
C3_l=dv_l_1^2
C3_s=dv_s_1^2
%plot parameters
axis equal;
%name=char(strcat('Transfer Opportunity',{' '},num2str(opp)));
%title(name)
xlabel('X distance (km)')
ylabel('Y distance (km)')
zlabel('Z distance (km)')
xlim([-1.5e+009 1.5e+009]);
ylim([-1.5e+009 1.5e+009]);
az = 0;
el = 90;
view(az, el);
legend([{'Sun'},{'Mercury'},{'Venus'},{'Earth'},{'Mars'},{'Jupiter'},{'Saturn'},{'Earth Departure'},{'Saturn Arrival'},{'Lambert Trajectory (long way)'},{'Lambert Trajectory (short way)'}]);

%% Porkchop Plot
% 
% tof=2*365; %days
% ToF=86400*tof; %seconds
% d_dep=[1:30:378];
% d_arr=d_dep+800;
% for i=d_dep
%     for j=d_arr
%         tof=j-i;
%         ToF=86400*tof;
%         r1=system_N(i,19:21);
%         r2=system_N(tof,37:39);
%         [v1,v2] = Lambert_solver(r1,r2,ToF,'long',mu(1));
%         dv_1=norm(system_N(i,22:24)-v1);
%         dv_2=norm(system_N(tof,40:42)-v2);
%         C3(i,j)=dv_1^2;
%     end
% end

%contour(d_dep,d_arr,C3)