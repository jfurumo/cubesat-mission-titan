%John Furumo
%ENAE 404 section 0101

%Patched Conics Worksheet

clear;
clf;
clc;
close all;

format shorteng;

mu_sun=1.327e+011;   %km^3/s^2

mu_earth=3.986e+005; %km^3/s^2
r_earth=149.5e+006;  %km

mu_mars=4.305e+004;  %km^3/s^2
r_mars=227.8e+006;   %km

mu_mercury=2.232e+004; %km^3/s^2
r_mercury=57.9e+006    %km

mu_jupiter=1.268e+008; %km^3/s^2
r_jupiter=778e+006;    %km

mu_s=mu_sun;
mu_e=mu_earth;
%mu_m=mu_mars;
mu_m=mu_mercury;
mu_j=mu_jupiter;

r_e=r_earth;
%r_m=r_mars;
r_m=r_mercury;
r_j=r_jupiter;

%Earth orbital velocity (circular orbit assumption)
v_earth=sqrt(mu_s/r_e) %km/s
%Spacecraft transfer ellipse velocity at periapsis
v_t_p=sqrt(((2*mu_s)/r_e)-(2*mu_s)/(r_e+r_m)) %km/s
%difference between spacecraft and earth velocity
v_inf_e=v_t_p-v_earth %km/s
%Mars orbital velocity (circular orbit assumption)
v_mars=sqrt(mu_s/r_m) %km/s
%Spacecraft transfer ellipse velocity at apoapsis
v_t_a=sqrt(((2*mu_s)/r_m)-(2*mu_s)/(r_e+r_m)) %km/s
%difference between mars and spacecraft velocity
v_inf_m=v_mars-v_t_a %km/s

%circular Earth parking orbit
r1=7000; %km
v1=sqrt(mu_e/r1) %km/s
%hyperbolic Earth escape orbit
v_h_1=sqrt(v_inf_e^2+(2*mu_e)/r1) %km/s
%Earth departure delta V
delta_v_1=v_h_1-v1 %km/s

%circular Mars parking orbit
r2=5000; %km
v2=sqrt(mu_m/r2) %km/s
%hyperbolic Mars capture orbit
v_h_2=sqrt(v_inf_m^2+(2*mu_m)/r2) %km/s
%Mars capture delta V
delta_v_2=v_h_2-v2 %km/s

%total delta V
delta_v_total=delta_v_1+delta_v_2 %km/s

