% Differential Equations for the bottle rocket

function datafunc = ode45func2(t,conditions_0)

%% Global Variables
global g discharge_coeff rho_air_amb volume_bottle p_amb gamma rho_water ...
    area_throat area_bottle gas_const_air cd x_0 z_0 l_stand mass_air_0 ...
    heading_0 volume_air_0 p_abs_0 windVelG windVelA;

%% Initial Conditions
x = conditions_0(1); % x position
y = conditions_0(2); % y position
z = conditions_0(3); % z position
vel_x = conditions_0(4); % Velocity downrange direction
vel_y = conditions_0(5); % Velocity crossrange direction
vel_z = conditions_0(6); % Velocity vertical direction
mass_rocket = conditions_0(7); % Mass of the bottle rocket
mass_air = conditions_0(8); % mass of air in the bottle rocket
volume_air = conditions_0(9); % volume of air inside the bottle rocket

%% Calculations
windVel = (windVelG + windVelA) / 2;
if sqrt((x - x_0)^(2) + (z - z_0)^(2)) > l_stand
        vel_rel = [vel_x-windVel(1), vel_y-windVel(2), vel_z-windVel(3)];
        vel_rel_tot = sqrt((vel_x + windVel(1))^2 + (vel_y + windVel(2))^2 + (vel_z + windVel(3))^2);
else
    vel_rel = [vel_x, vel_y, vel_z];
    vel_rel_tot = sqrt((vel_x)^2 + (vel_y)^2 + (vel_z)^2);
end
    
drag = (rho_air_amb / 2) * (vel_rel_tot)^2 * cd * area_bottle; % find drag using equation 2 [N]

    if volume_air >= volume_bottle % if the volume of the air is greater than or equal to the volume of the bottle (signifies the end of the water thrust stage)

        p_water_thrust_end = p_abs_0 * ...
            ((volume_air_0 / volume_bottle)^(gamma)); % pressure at the end of the water thrust stage, using equation 13

        p_air_thrust_stage = p_water_thrust_end * ...
            ((mass_air / mass_air_0)^(gamma)); % pressure of air during air thrust stage, using equation 14

        rho_air_thrust_stage = mass_air / volume_bottle; % denstiy of air at p_air_thrust_stage, using equation 15
        temp_air_thrust_stage = p_air_thrust_stage / ...
            (rho_air_thrust_stage * gas_const_air); % temperature of air at p_air_thrust_stage, using equation 15

    end

    if volume_air < volume_bottle % if the volume of air is less than the volume of the bottle (signifies the bottle rocket is still in the water thrust stage)

        p_water_thrust_stage = p_abs_0 * (volume_air_0 / volume_air)^(gamma); % pressure in the bottle, using equation 3

        thrust = 2 * discharge_coeff * area_throat * ...
            (p_water_thrust_stage - p_amb); % thrust of the bottle rocket [N], using equation 5

        dvolumeair_dt = discharge_coeff * area_throat * sqrt((2 / rho_water) * ...
            (p_abs_0 * ((volume_air_0 / volume_air)^(gamma)) - p_amb)); % change in volume over change in time, using equation 9

        dmasstot_dt = - discharge_coeff * rho_water * area_throat * ...
            sqrt((2 * (p_water_thrust_stage - p_amb)) / rho_water); % change in the total mass of the bottle rocket over change in time, using equation 10

        dmassair_dt = 0; % change in mass of air over change in time, using equation 24 (zero for the water thrust stage)

    elseif p_air_thrust_stage > p_amb % if p_air_thrust_stage is greater than p_amb (signifies the bottle rocket is in the air thrust phase)

        p_critical = p_air_thrust_stage * ...
            ((2 / (gamma + 1))^(gamma / (gamma - 1))); % critical pressure in the bottle, using equation 16

        if p_critical > p_amb % if p_critical is greater than p_amb (checking for chocked flow)

          mach_exit = 1; % definition of chocked flow

          temp_exit = (2 / (gamma + 1)) * temp_air_thrust_stage; % exit temperature, using equation 18
          p_exit = p_critical; % exit pressure equals critical pressure if at chocked flow, using equation 18
          rho_exit = p_exit / (gas_const_air * temp_exit); % exit density, using equation 18

          vel_exit = sqrt(gamma * gas_const_air * temp_exit); % exit velocity, using equation 17

        else % if not chocked flow

            mach_exit = sqrt((2 / (gamma - 1)) * ...
                (((p_air_thrust_stage / p_amb)^((gamma - 1) / gamma)) - 1)); % exit mach number, using equation 19
            temp_exit = temp_air_thrust_stage / (1 + ((gamma - 1) / 2) * ...
                (mach_exit)^(2)); % exit temperature, using equation 20
            p_exit = p_amb; % exit pressure equals ambient pressure, using equation 20
            rho_exit = p_exit / (gas_const_air * temp_exit); % exit density, using equation 20

            vel_exit = mach_exit * sqrt(gamma * gas_const_air * temp_exit); % exit velocity, using equation 21

        end   
        
        dmassair_dt = -1 * ...
            (discharge_coeff * rho_exit * area_throat * vel_exit); % change in mass of the air in the rocket over change in time, using equation 23
        
        thrust = (-1 * dmassair_dt * vel_exit) + (p_amb - p_exit) * ...
            area_throat; % thrust of the rocket, using equation 22
        
        dmasstot_dt = dmassair_dt; % change in total mass of the bottle rocket over time, using equation 24 
        
        dvolumeair_dt = 0; % change in volme of air over change in time (zero because the voulme of air no longer changes since the bottle rocket is in air thrust stage)
        
    else % if the bottle rocket is in stage 3 flight, no thrust
        
        thrust = 0; % no thrust in stage 3
        dvolumeair_dt = 0; % volume of air is not changing in stage 3
        dmasstot_dt = 0; % mass of bottle rocket is not changing in stage 3
        dmassair_dt = 0; % mass of air is not changing in stage 3
        
    end

    % x, y, and z components of velocity
    dx_dt = vel_x;
    dy_dt = vel_y;
    dz_dt = vel_z;
    
    if sqrt((x - x_0)^(2) + (z - z_0)^(2)) < l_stand % if the bottle rocket hasn't lef the test stand  
        
        % initial headings, component wise 
        x_height = heading_0(1);
        y_height = heading_0(2);
        z_height = heading_0(3);
        
    else % if the bottle rocket has left the test stand
        
        % heading, component wise
        x_height = vel_rel(1) / abs(vel_rel_tot);
        y_height = vel_rel(2) / abs(vel_rel_tot);
        z_height = vel_rel(3) / abs(vel_rel_tot);
        
    end
    
    % acceleration in each component
    dvelx_dt = (thrust - drag) * (x_height / mass_rocket);
    dvely_dt = (thrust - drag) * (y_height / mass_rocket);
    dvelz_dt = (thrust - drag) * (z_height / mass_rocket) - g;
  
    % differentials for ode45
    datafunc(1) = dx_dt;
    datafunc(2) = dy_dt;
    datafunc(3) = dz_dt;
    datafunc(4) = dvelx_dt;
    datafunc(5) = dvely_dt;
    datafunc(6) = dvelz_dt;
    datafunc(7) = dmasstot_dt;
    datafunc(8) = dmassair_dt;
    datafunc(9) = dvolumeair_dt;
    datafunc = datafunc'; % get equationResults into a column vec for ode45 
    
end