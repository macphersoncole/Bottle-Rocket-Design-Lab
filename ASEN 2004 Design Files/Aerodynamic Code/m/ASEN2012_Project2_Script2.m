%% ASEN 2012 Project 2 (Water Bottle Rocket)

% Authors
%
% 1) Cole MacPherson
% 2) Ankrit Uprety
%
% https://www.timeanddate.com/weather/usa/boulder/historic?month=3&year=2020

%% Housekeeping

clc;
clear;
close all;

%% Declare Global Variables and Constants/Verification Constants

% Gloabl Variables Declaration

global g discharge_coeff rho_air_amb volume_bottle p_amb gamma rho_water ...
    diameter_throat area_throat area_bottle diameter_bottle gas_const_air ...
    mass_bottle cd p_gage_0 volume_water_0 temp_air_0 vel_0 theta_0 x_0 y_0 ...
    z_0 vel0_x vel0_y vel0_z l_stand mass_air_0 mass_rocket_0 heading_0 ...
    volume_air_0 p_abs_0 tspan windVelG windVelA;

% Constants

g = 9.81; % acceleration due to gravity [m/s^2]
discharge_coeff = 0.8; % discharge coefficient
rho_air_amb = 1.17358; % ambient air density [kg/m^3] GOLD ROCKET -> 1.17358 [kg/m^3] TA BASELINE -> 1.23478
volume_bottle = 0.002; % volume of the empty bottle [m^3]
p_amb = 99640; % atmospheric pressure [pa] GOLD ROCKET -> 99640 [pa] TA BASELINE -> 101592
gamma = 1.4; % ratio of specific heats for heat
rho_water = 1000; % density of water [kg/m^2]
diameter_throat = 2.1 / 100; % diameter of throat [m]
diameter_bottle = 10.5 / 100; % diameter of bottle [m]
area_throat = (1/4) * pi * (diameter_throat)^2; % area of throat [m^2]
area_bottle = (1/4) * pi * (diameter_bottle)^2; % area of bottle [m^2]
gas_const_air = 287; % gas constant of air [J/(kg*K)]
mass_bottle = 0.126; % mass of the empty bottle w/ cone and fins [kg]
cd = 0.327; % coefficient of drag
l_stand = 0.5; % langth of test stand [m]
tspan = [0, 5]; % integration time span [s]

% Max Distance Constants

p_gage_0 = 40.5 * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
volume_water_0 = 0.0005833; % intitial volume of water inside bottle [m^3] 
volume_air_0 = volume_bottle - volume_water_0; % intitial volume of air inside the bottle [m^3]
temp_air_0 = 295.15; % initial temperature of air [K] GOLD ROCKET -> 295.15 TA BASELINE -> 283.15
vel_0 = 0; % initial velocity of bottle rocket [m/s]
vel0_x = 0; % initial velocity of bottle rocket in the x direction [m/s]
vel0_y = 0; % initial velocity of bottle rocket in the y direction [m/s]
vel0_z = 0; % initial velocity of bottle rocket in the z direction [m/s]
theta_0 = 45; % intitial angle of bottle rocket [degrees] MAXIMUM DISTANCE -> 40
x_0 = 0; % initial downrange distance [m]
y_0 = 0; % initial crossrange distance [m]
z_0 = 0.25; % initial verital height [m]
windTheta = 0;
windSPDG = 0; % wind velocity [m/s]
windSPDA = 0;
windDir = [cosd(windTheta) sind(windTheta) 0]; % wind velocity direction [m/s]
windVelG = [windDir(1)*windSPDG, windDir(2)*windSPDG, windDir(3)*windSPDG];
windVelA = [windDir(1)*windSPDA, windDir(2)*windSPDA, windDir(3)*windSPDA];
heading_0 = [cosd(theta_0),0, sind(theta_0)]; % initial heading of bottle rocket [x,y]
mass_air_0 = (p_abs_0 * volume_air_0) / (gas_const_air * temp_air_0); % initial mass of air [kg]
mass_rocket_0 = mass_bottle + (rho_water * volume_water_0) + mass_air_0; % total mass of the bottle rocket [kg]


%% Calculations

initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

[time, data] = ode45('ode45func2',tspan,initial_conditions);

%[thrust, phase1_end_i, phase2_end_i] =  thrustvecfunc2(time,data);

%% Plot trajectory of bottle rocket
i = 1;
while data(i,3) >= 0
    i = i+1;
end

max_distance_x = (data(i-1,1)+data(i,1)) / 2;
max_distance_y = (data(i-1,2)+data(i,2)) / 2;
max_distance_z = 0;
max_distances = [max_distance_x,max_distance_y,max_distance_z];
distanceFromLaunch = sqrt(max_distance_x^2+max_distance_y^2);

x_limit = max_distances(1)*1.05;
y_limit = max_distances(2)*1.05;
z_limit = max(data(:,3))*1.05;

j = 1;
while data(j,3) ~= max(data(:,3))
    j = j+1;
    max_height = [data(j,1) data(j,2) data(j,3)];
end

figure(1)
plot3(x_0,y_0,z_0,'or',max_distances(1),max_distances(2), ...
    max_distances(3),'*r',max_height(1),max_height(2),max_height(3),'+r', ...
    data(1:i,1),data(1:i,2),data(1:i,3),'b','LineWidth',2);
xlim([x_0 x_limit]);
if y_limit > 0
    ylim([y_0 y_limit]);
elseif y_limit < 0
    ylim([y_limit y_0]);
end
zlim([z_0-0.3 z_limit]);
title('Rocket Trajectory');
xlabel('Downrange Distance, [m]');
ylabel('Crossrange Distance, [m]');
zlabel('Vertical Distance, [m]');
legend(['Launch Stand Location'],['Maximum Downrange Distance, ' num2str(max_distances(1)) ' [m]'],['Maximum Height, ' num2str(max_height(3)) ' [m]'],'location','northoutside');
grid on;

%% Montecarlo Error Ellipses

xyCoorSim = zeros(100,2);

for i = 1:101
    
    rho_air_amb = 1.17358 + (0.1*rand - 0.05); % ambient air density [kg/m^3] GOLD ROCKET -> 1.17358 [kg/m^3] TA BASELINE -> 1.23478
    p_amb = 99640 + (10*rand - 5); % atmospheric pressure [pa] GOLD ROCKET -> 99640 [pa] TA BASELINE -> 101592
    cd = 0.327 + 0.01625*randn; % coefficient of drag
    diameter_throat = (2.1 / 100) + 0.0005*randn; % diameter of throat [m]
    diameter_bottle = (10.5 / 100) + 0.0005*randn; % diameter of bottle [m]
    area_throat = (1/4) * pi * (diameter_throat)^2; % area of throat [m^2]
    area_bottle = (1/4) * pi * (diameter_bottle)^2; % area of bottle [m^2]
    
    p_gage_0 = (40.5 + 0.5*randn) * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
    p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
    volume_water_0 = 0.0005833 + (0.0000005*randn); % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
    volume_air_0 = volume_bottle - volume_water_0; % intitial volume of air inside the bottle [m^3]
    temp_air_0 = 295.15 + (1*rand - 0.5); % initial temperature of air [K] GOLD ROCKET -> 295.15 TA BASELINE -> 283.15
    theta_0 = 45 + randn; % intitial angle of bottle rocket [degrees] MAXIMUM DISTANCE -> 40
    windTheta = 360*rand;
    windSPDG = (10*rand)/2.237; % wind ground velocity [m/s]
    windSPDA = (10*rand)/2.237; % wind aloft velocity [m/s] 
    windDir = [cosd(windTheta) sind(windTheta) 0]; % wind velocity direction [m/s]
    windVelG = [windDir(1)*windSPDG, windDir(2)*windSPDG, windDir(3)*windSPDG];
    windVelA = [windDir(1)*windSPDA, windDir(2)*windSPDA, windDir(3)*windSPDA];
    heading_0 = [cosd(theta_0),0, sind(theta_0)]; % initial heading of bottle rocket [x,y]
    mass_air_0 = (p_abs_0 * volume_air_0) / (gas_const_air * temp_air_0); % initial mass of air [kg]
    mass_rocket_0 = mass_bottle + (rho_water * volume_water_0) + mass_air_0; % total mass of the bottle rocket [kg]

    
    initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

    [timeSim, dataSim] = ode45('ode45func2',tspan,initial_conditions);
    
    j = 1;
    while dataSim(j,3) >= 0
        j = j+1;
    end
    j = j-1;
    
    xyCoorSim(i,:) = [(dataSim(j,1)+dataSim(j+1,1))/2, (dataSim(j,2)+dataSim(j+1,2))/2];
    
end

error_ellipses(xyCoorSim,0);

%% Montecarlo Error Ellipses GOLD ROCKET

xyCoorSim_g = zeros(100,2);

for i = 1:101
    
    rho_air_amb = 1.17358 + (0.1*rand - 0.05); % ambient air density [kg/m^3] GOLD ROCKET -> 1.17358 [kg/m^3] TA BASELINE -> 1.23478
    p_amb = 99640 + (1000*rand - 500); % atmospheric pressure [pa] GOLD ROCKET -> 99640 [pa] TA BASELINE -> 101592
    cd = 0.327 + 0.004*randn; % coefficient of drag
    diameter_throat = (2.1 / 100) + 0.0005*randn; % diameter of throat [m]
    diameter_bottle = (10.5 / 100) + 0.0005*randn; % diameter of bottle [m]
    area_throat = (1/4) * pi * (diameter_throat)^2; % area of throat [m^2]
    area_bottle = (1/4) * pi * (diameter_bottle)^2; % area of bottle [m^2]
    
    p_gage_0 = (40.5 + (1*rand - 0.5)) * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
    p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
    volume_water_0 = 0.000993 + (0.0000005*randn); % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
    volume_air_0 = volume_bottle - volume_water_0; % intitial volume of air inside the bottle [m^3]
    temp_air_0 = 295.15 + (2*rand - 1); % initial temperature of air [K] GOLD ROCKET -> 295.15 TA BASELINE -> 283.15
    theta_0 = 45 + randn; % intitial angle of bottle rocket [degrees] MAXIMUM DISTANCE -> 40
    windTheta = (136+(22.5*rand-11.25)) + 1.5*randn;
    windSPDG = 0 + (2*0.44704*rand-0.44704); % wind ground velocity [m/s]
    windSPDA = 0.44704 + (2*0.44704*rand-0.44704); % wind aloft velocity [m/s] 
    windDir = [cosd(windTheta) sind(windTheta) 0]; % wind velocity direction [m/s]
    windVelG = [windDir(1)*windSPDG, windDir(2)*windSPDG, windDir(3)*windSPDG];
    windVelA = [windDir(1)*windSPDA, windDir(2)*windSPDA, windDir(3)*windSPDA];
    heading_0 = [cosd(theta_0),0, sind(theta_0)]; % initial heading of bottle rocket [x,y]
    mass_air_0 = (p_abs_0 * volume_air_0) / (gas_const_air * temp_air_0); % initial mass of air [kg]
    mass_rocket_0 = mass_bottle + (rho_water * volume_water_0) + mass_air_0; % total mass of the bottle rocket [kg]

    
    initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

    [timeSim_g, dataSim_g] = ode45('ode45func2',tspan,initial_conditions);
    
    j = 1;
    while dataSim_g(j,3) >= 0
        j = j+1;
    end
    j = j-1;
    
    xyCoorSim_g(i,:) = [(dataSim_g(j,1)+dataSim_g(j+1,1))/2, (dataSim_g(j,2)+dataSim_g(j+1,2))/2];
    
end

error_ellipses(xyCoorSim_g,1);


    


