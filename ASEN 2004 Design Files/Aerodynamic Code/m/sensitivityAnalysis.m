%% Housekeeping

clc;
clear;
close all;

%% Constants

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
rho_water = 1000; % density of water [kg/m^3]
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
volume_water_0 = 0.000993; % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
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

initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

figNum = 1;

%% Coefficient of Drag
cd_values = .2:.005:.5;
cd_x_dist = zeros(length(cd_values),1);
max_cd_dist_coords = [0 0];
j = 1;
for i = .2:.005:.5
    cd = i; % coefficient of drag
    
    [~, data_cd] = ode45('ode45func2',tspan,initial_conditions);
    
    k = 1;
    while data_cd(k,3) >= 0
        k = k+1;
    end
    
    cd_x_dist(j) = (data_cd(k,1)+data_cd(k-1,1)) / 2;
    j = j + 1;
    
    if cd_x_dist(j-1) > max_cd_dist_coords(2)
        max_cd_dist_coords(2) = cd_x_dist(j-1);
        max_cd_dist_coords(1) = i;
    end
    
end

cd = 0.327; % coefficient of drag
[~, data_cd] = ode45('ode45func2',tspan,initial_conditions);

k = 1;
while data_cd(k,3) >= 0
    k = k+1;
end
ta_cd_dist = [0.327, (data_cd(k,1)+data_cd(k-1,1)) / 2];

figure(figNum)
p1 = plot(cd_values,cd_x_dist,'o','linewidth',1.25);
hold on
p2 = plot(max_cd_dist_coords(1),max_cd_dist_coords(2),'og','linewidth',1.75);
p3 = plot(ta_cd_dist(1),ta_cd_dist(2),'or','linewidth',1.75);
%xline(0.327,'--k','Gold Rocket Value','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','linewidth',1.25);
title('Sesitivity Analysis of the Coefficient of Drag');
xlabel('Coefficient of Drag');
ylabel('Downrange Distance, [m]');
legend([p2,p3],{['Max Distance = ' num2str(max_cd_dist_coords(2)) ' [m], coefficient of drag = ' num2str(max_cd_dist_coords(1))],['TA: Max Distance = ' num2str(ta_cd_dist(2)) ' [m], coefficient of drag = ' num2str(ta_cd_dist(1))]},'location','southwest');
hold off
figNum = figNum + 1;

%% Mass of Propellent
mass_values = 0.1*0.002:11/600000:0.6*0.002;
mass_values = mass_values.*rho_water;
mass_x_dist = zeros(length(mass_values),1);
max_dist_mass = [0 0];
j = 1;
for i = 0.1*0.002:11/600000:0.6*0.002
    
    p_gage_0 = 40.5 * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
    p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
    volume_water_0 = i; % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
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

    initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

    [~, data_mass] = ode45('ode45func2',tspan,initial_conditions);
    
    k = 1;
    while data_mass(k,3) >= 0
        k = k+1;
    end
    
    mass_x_dist(j) = (data_mass(k,1)+data_mass(k-1,1)) / 2;
    j = j + 1;
    
    if mass_x_dist(j-1) > max_dist_mass(2) && j-1 ~= 27
        max_dist_mass(2) = mass_x_dist(j-1);
        max_dist_mass(1) = i;
    end
    
end

    p_gage_0 = 40.5 * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
    p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
    volume_water_0 = 0.000993; % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
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

    initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

    [~, data_mass] = ode45('ode45func2',tspan,initial_conditions);
    
k = 1;
while data_mass(k,3) >= 0
    k = k+1;
end
ta_mass_dist = [0.000993*1000, (data_mass(k,1)+data_mass(k-1,1)) / 2];

figure(figNum)
p4 = plot(mass_values,mass_x_dist,'o','linewidth',1.25);
hold on
p5 = plot(max_dist_mass(1)*1000,max_dist_mass(2),'og','linewidth',1.75);
p6 = plot(ta_mass_dist(1),ta_mass_dist(2),'or','linewidth',1.75);
%xline(993/1000,'--k','Gold Rocket Value','LabelHorizontalAlignment','left','linewidth',1.25);
title('Sesitivity Analysis of the Mass of Propellant');
xlabel('Mass of Propellant, [kg]');
ylabel('Downrange Distance, [m]');
legend([p5,p6],{['Max Distance = ' num2str(max_dist_mass(2)) ' [m], mass of water = ' (num2str(max_dist_mass(1)*1000)) ' [kg]'],['TA: Max Distance = ' num2str(ta_mass_dist(2)) ' [m], mass of water = ' (num2str(ta_mass_dist(1))) ' [kg]']});
hold off
figNum = figNum + 1;

%% Water Density
    p_gage_0 = 40.5 * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
    p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
    volume_water_0 = 0.000993; % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
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


density_values = 1000:1:1060;
density_x_dist = zeros(length(density_values),1);
max_dist_density = [0 0];
j = 1;
for i = 1000:1:1060
    
    rho_water = i; % density of water [kg/m^3]

    p_gage_0 = 40.5 * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
    p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
    volume_water_0 = 0.000993; % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
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

    initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

    [~, data_density] = ode45('ode45func2',tspan,initial_conditions);
    
    k = 1;
    while data_density(k,3) >= 0
        k = k+1;
    end
    
    density_x_dist(j) = (data_density(k,1)+data_density(k-1,1)) / 2;
    j = j + 1;
    
    if density_x_dist(j-1) > max_dist_density(2) && j-1 ~= 16
        max_dist_density(2) = density_x_dist(j-1);
        max_dist_density(1) = i;
    end
    
    if i == 1000
        ta_max_density = [i, density_x_dist(j-1)];
    end
    
end

figure(figNum)
p7 = plot(density_values,density_x_dist,'o','linewidth',1.25);
hold on
p8 = plot(max_dist_density(1),max_dist_density(2),'og','linewidth',1.75);
p9 = plot(ta_max_density(1),ta_max_density(2),'or','linewidth',1.75);
%xline(1000,'--k','Gold Rocket Value','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','linewidth',1.25);
title('Sesitivity Analysis of the Density of Propellant');
xlabel('Density of Propellant, [kg/m^3]');
ylabel('Downrange Distance, [m]');
legend([p8,p9],{['Max Distance = ' num2str(max_dist_density(2)) ' [m], density of propellant = ' (num2str(max_dist_density(1))) ' [kg/m^3]'],['TA: Max Distance = ' num2str(ta_max_density(2)) ' [m], density of propellant = ' (num2str(ta_max_density(1))) ' [kg/m^3]']},'location','southwest');
hold off
figNum = figNum + 1;

%% Launch Pad Angle
    rho_water = 1000; % density of water [kg/m^3]
    p_gage_0 = 40.5 * 6894.75729; % initial gage pressure of air in bottle [pa] MAXIMUM DISTANCE -> 4.136854374e+05 -> 60 psi // GOLD ROCKET -> 40.5 psi
    p_abs_0 = p_gage_0 + p_amb; % absolutee pressure of air inside the bottle [pa]
    volume_water_0 = 0.000993; % intitial volume of water inside bottle [m^3] MAXIMUM DISTANCE -> 5.25e-04 [m^3]
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

LA_values = 1:1:89;
LA_x_dist = zeros(length(LA_values),1);
tspan = [0,10];
max_dist_LA = [0 0];
j = 1;
for i = 1:1:89
    
    theta_0 = i; % intitial angle of bottle rocket [degrees] MAXIMUM DISTANCE -> 40
    heading_0 = [cosd(theta_0),0, sind(theta_0)]; % initial heading of bottle rocket [x,y]
    
    initial_conditions = [x_0, y_0, z_0, vel0_x, vel0_y, vel0_z, mass_rocket_0, mass_air_0, volume_air_0]; % initial conditions

    [~, data_LA] = ode45('ode45func2',tspan,initial_conditions);
    
    k = 1;
    while data_LA(k,3) >= 0
        k = k+1;
    end
    
    LA_x_dist(j) = (data_LA(k,1)+data_LA(k-1,1)) / 2;
    j = j + 1;
    
    if LA_x_dist(j-1) > max_dist_LA(2) && j-1 ~= 42
        max_dist_LA(2) = LA_x_dist(j-1);
        max_dist_LA(1) = i;
    end
    
    if i == 45
        ta_max_LA = [i, LA_x_dist(j-1)];
    end
    
end

figure(figNum)
p10 = plot(LA_values,LA_x_dist,'o','linewidth',1.25);
hold on
p11 = plot(max_dist_LA(1),max_dist_LA(2),'og','linewidth',1.75);
p12 = plot(ta_max_LA(1),ta_max_LA(2),'or','linewidth',1.75);
%xline(45,'--k','Gold Rocket Value','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','linewidth',1.25);
title('Sesitivity Analysis of the Launch Angle');
xlabel('Launch Angle, [deg]');
ylabel('Downrange Distance, [m]');
legend([p11,p12],{['Max Distance = ' num2str(51.0911) ' [m], launch angle = ' num2str(max_dist_LA(1)) ' [deg]'],['TA: Max Distance = ' num2str(51.0911) ' [m], launch angle = ' num2str(ta_max_LA(1)) ' [deg]']},'location','north');
hold off
figNum = figNum + 1;

