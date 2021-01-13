%% ASEN 2012 Project 2 (Water Bottle Rocket)

% Authors
%
% 1) Cole MacPherson
% 2) Ankrit Uprety
%

%% Housekeeping

clc;
clear;
close all;

%% Declare Global Variables and Constants/Verification Constants

% Gloabl Variables Declaration

global g discharge_coeff rho_air_amb volume_bottle p_amb gamma rho_water ...
    diameter_throat area_throat area_bottle diameter_bottle gas_const_air ...
    mass_bottle cd p_gage_0 volume_water_0 temp_air_0 vel_0 theta_0 x_0 y_0 ...
    l_stand mass_air_0 mass_rocket_0 heading_0 volume_air_0 p_abs_0 tspan;

% Constants

g = 9.81; % acceleration due to gravity [m/s^2]
discharge_coeff = 0.8; % discharge coefficient
rho_air_amb = 0.961; % ambient air density [kg/m^3]
volume_bottle = 0.002; % volume of the empty bottle [m^3]
p_amb = 12.1 * 6894.75729; % atmospheric pressure [pa]
gamma = 1.4; % ratio of specific heats for heat
rho_water = 1000; % density of water [kg/m^2]
diameter_throat = 2.1 / 100; % diameter of throat [m]
diameter_bottle = 10.5 / 100; % diameter of bottle [m]
area_throat = (1/4) * pi * (diameter_throat)^2; % area of throat [m^2]
area_bottle = (1/4) * pi * (diameter_bottle)^2; % area of bottle [m^2]
gas_const_air = 287; % gas constant of air [J/(kg*K)]
mass_bottle = 0.15; % mass of the empty bottle w/ cone and fins [kg]
cd = 0.5; % coefficient of drag
l_stand = 0.5; % langth of test stand [m]
tspan = [0, 5]; % integration time span [s]

% Verification Constants
%{
p_gage_0 = 50 * 6894.75729; % initial gage pressure of air in bottle [pa]
p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
volume_water_0 = 0.001; % intitial volume of water inside bottle [m^3]
volume_air_0 = volume_bottle - volume_water_0; % intitial volume of air inside the bottle [m^3]
temp_air_0 = 300; % intitial temperature of air [K]
vel_0 = 0; % intitial velocity of bottle rocket [m/s]
vel_x = 0; % intitial velocity of bottle rocket in the x direction [m/s]
vel_y = 0; % intitial velocity of bottle rocket in the y direction [m/s]
theta_0 = 45; % intitial angle of bottle rocket [degrees]
x_0 = 0; % intitial horizontal distance [m]
y_0 = 0.25; % initial vertical height [m]
heading_0 = [cosd(theta_0), sind(theta_0)]; % initial heading of bottle rocket [x,y]
mass_air_0 = (p_abs_0 * volume_air_0) / (gas_const_air * temp_air_0); % initial mass of air [kg]
mass_rocket_0 = mass_bottle + (rho_water * volume_water_0) + mass_air_0; % total mass of the bottle rocket [kg]
%}

% 80 meter Constants

p_gage_0 = 62.5 * 6894.75729; % initial gage pressure of air in bottle [pa]
p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
volume_water_0 = 0.0005; % intitial volume of water inside bottle [m^3]
volume_air_0 = volume_bottle - volume_water_0; % intitial volume of air inside the bottle [m^3]
temp_air_0 = 300; % intitial temperature of air [K]
vel_0 = 0; % intitial velocity of bottle rocket [m/s]
vel_x = 0; % intitial velocity of bottle rocket in the x direction [m/s]
vel_y = 0; % intitial velocity of bottle rocket in the y direction [m/s]
theta_0 = 40; % intitial angle of bottle rocket [degrees]
x_0 = 0; % intitial horizontal distance [m]
y_0 = 0.25; % initial vertical height [m]
heading_0 = [cosd(theta_0), sind(theta_0)]; % initial heading of bottle rocket [x,y]
mass_air_0 = (p_abs_0 * volume_air_0) / (gas_const_air * temp_air_0); % initial mass of air [kg]
mass_rocket_0 = mass_bottle + (rho_water * volume_water_0) + mass_air_0; % total mass of the bottle rocket [kg]

% initializing end times to assist in plotting
phase1_end_i = 0;
phase2_end_i = 0;

%% Calculations

initial_conditions = [x_0, y_0, vel_x, vel_y, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

[time, data] = ode45('ode45func',tspan,initial_conditions);

[thrust, phase1_end_i, phase2_end_i] =  thrustvecfunc(time,data);

% plot trajectory of bottle rocket
figure(1)
plot(data(:,1),data(:,2));
xline(80,'--',{'Target Distance'});
%{
xline(data(phase1_end_i,1),'--',{'Phase 1 Ending'});
xline(data(phase2_end_i,1),'--',{'Phase 2 Ending'});
%}
xlim([0 85]);
ylim([0 25]);
title('Rocket Trajectory');
xlabel('Horizontal Distance, [m]');
ylabel('Vertical Distance, [m]');

% plot x and y velocity of bottle rocket
figure(2)
plot(data(:,3),data(:,4));
xline(data(phase1_end_i,3),'--');
xline(data(phase2_end_i,3),'--');
title('Rocket Velocity');
xlabel('Horizontal Velocity, [m/s]');
ylabel('Vertical Velocity, [m/s]');

% plot x velocity vs time
figure(3)
plot(time,data(:,3));
xline(time(phase1_end_i,1),'--');
xline(time(phase2_end_i,1),'--');
title('Rocket Horizontal Velocity');
xlabel('Time, [s]');
ylabel('Horizontal Velocity, [m/s]');

% plot y velocity vs time
figure(4)
plot(time,data(:,4));
xline(time(phase1_end_i,1),'--');
xline(time(phase2_end_i,1),'--');
title('Rocket Vertical Velocity');
xlabel('Time, [s]');
ylabel('Vertical Velocity, [m/s]');

% plot total mass vs time
figure(5)
plot(time,data(:,5));
xline(time(phase1_end_i,1),'--',{'Phase 1 Ending'});
xline(time(phase2_end_i,1),'--',{'Phase 2 Ending'});
title('Rocket Total mass');
xlim([0 0.25]);
xlabel('Time, [s]');
ylabel('Total Mass, [kg]');

% plot mass of air vs time
figure(6)
plot(time,data(:,6));
xline(time(phase1_end_i,1),'--');
xline(time(phase2_end_i,1),'--');
xlim([0 0.05]); 
title('Mass of Air in Rocket vs Time');
xlabel('Time, [s]');
ylabel('Mass of Air in Rocket, [kg]');

% plot volume of air vs time
figure(7)
plot(time,data(:,7));
xline(time(phase1_end_i,1),'--',{'Phase 1','Ending'});
xline(time(phase2_end_i,1),'--',{'Phase 2','Ending'});
xlim([0 .0825]); 
title('Volume of Air');
xlabel('Time, [s]');
ylabel('Vomume of Air, [m^3]');



% plot thrust vs time
figure(8)
plot(time,thrust); 
xline(time(phase1_end_i,1),'--',{'Phase 1 Ending'});
xline(time(phase2_end_i,1),'--',{'Phase 2 Ending'});
xlim([0 .45]); 
title('Thrust');
xlabel('Time, [s]');
ylabel('Thrust, [N]');
