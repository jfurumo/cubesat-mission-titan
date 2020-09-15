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
G=6.67428E-020; %km^3/kg*s^2
Spacecraft_mass=100; %4x 12U cubesats of 25kg each
n=9; %number of bodies: Saturn + 8 moons

%Saturnian system body
body{1}=('Saturn.txt');
body{2}=('Mimas.txt');
body{3}=('Enceladus.txt');
body{4}=('Tethys.txt');
body{5}=('Dione.txt');
body{6}=('Rhea.txt');
body{7}=('Titan.txt');
body{8}=('Hyperion.txt');
body{9}=('Iapetus.txt');

%Saturnian system mu values (km^3/s^2)
mu(1)=37931207.8; %km^3/s^2   %saturn     
mu(2)=2.504; %km^3/s^2        %mimas
mu(3)=7.211; %km^3/s^2        %enceladus
mu(4)=41.21; %km^3/s^2        %tethys
mu(5)=73.113; %km^3/s^2       %dione
mu(6)=153.94; %km^3/s^2       %rhea
mu(7)=8978.14; %km^3/s^2      %titan
mu(8)=0.3708; %km^3/s^2       %hyperion
mu(9)=120.53; %km^3/s^2       %iapetus
%switch to column vector
mu=mu';

for p=1:n
%folder pathway
folder=('C:\Users\jfuru\Desktop\UMD\Fall 2019\ENAE 601\final project\Solar System Ephemeris\2020-2034\Saturn Moons\Cronocentric\Barycenter\');
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
%end date:   January 1, 2026 (2192 days later)
G=6.67428E-020; %km^3/kg*s^2
%tf=n_ephemeris*86400; %s
tf=2*365*86400; %s
%one state vector per day
tspan=[0:86400:tf]; %s
%tspan=[0:1:n_ephemeris]; %day
tolerance=1e-006;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);


%Saturnian System initial state vector
for i=1:n
SS_state(1,(6*(i-1)+1):(6*i))=cell2mat((SS{i}(2)));
end
SS_state=SS_state';

%Numerical integration of solar system
[t,system_N]=ode45(@integrate_NBP,tspan,SS_state,options,n,mu);

%3D plot of solar system
figure('Name','Saturnian System');
plot3(system_N(:,1),system_N(:,2),system_N(:,3),'o','linewidth',2)   %saturn
hold on;
plot3(system_N(:,7),system_N(:,8),system_N(:,9),'linewidth',2)    %mimas
plot3(system_N(:,13),system_N(:,14),system_N(:,15),'linewidth',2) %enceladus
plot3(system_N(:,19),system_N(:,20),system_N(:,21),'linewidth',2) %tethys
plot3(system_N(:,25),system_N(:,26),system_N(:,27),'linewidth',2) %dione
plot3(system_N(:,31),system_N(:,32),system_N(:,33),'linewidth',2) %rhea
plot3(system_N(:,37),system_N(:,38),system_N(:,39),'linewidth',2) %titan
plot3(system_N(:,43),system_N(:,44),system_N(:,45),'linewidth',2) %hyperion
plot3(system_N(:,49),system_N(:,50),system_N(:,51),'linewidth',2) %iapetus
axis equal;
title({'Saturnian System','NBP Simulation 01/01/2034 - 01/01/2036'})
legend([{'Saturn'},{'Mimas'},{'Enceladus'},{'Tethys'},{'Dione'},{'Rhea'},{'Titan'},{'Hyperion'},{'Iapetus'}])
xlabel('X distance (km)')
ylabel('Y distance (km)')
zlabel('Z distance (km)')
% xlim([-1.5e+009 1.5e+009]);
% ylim([-.5e+009 1.5e+009]);
% az = 0;
% el = 0;
% view(az, el);
