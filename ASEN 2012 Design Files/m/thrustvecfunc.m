function [thrustvec, phase1_end, phase2_end] = thrustvecfunc(t,data)

    global discharge_coeff volume_bottle p_amb gamma area_throat ...
        gas_const_air mass_air_0 volume_air_0 p_abs_0;

    phase1_end_i = 0; % to plot line at stage endings
    phase2_end_i = 0; % to plot line at stage endings
    thrust = zeros(1,length(t)); % initializing the thrust vector
    
    for i = 1:length(t)
        
        mass_air = data(i,6); % mass of the air
        volume_air = data(i,7); % volume of air
        
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

            thrust(i) = 2 * discharge_coeff * area_throat * ...
                (p_water_thrust_stage - p_amb); % thrust of the bottle rocket [N], using equation 5

        elseif p_air_thrust_stage > p_amb % if p_air_thrust_stage is greater than p_amb (signifies the bottle rocket is in the air thrust phase)

            if phase1_end_i == 0
                
                phase1_end_i = i-1; % finds the indice for the end of phase 1
                
            end
            
            p_critical = p_air_thrust_stage * ...
                ((2 / (gamma + 1))^(gamma / (gamma - 1))); % critical pressure in the bottle, using equation 16

            if p_critical > p_amb % if p_critical is greater than p_amb (checking for chocked flow, definently won't happen but you asked for it)

              mach_exit = 1; % definition of chocked flow

              p_exit = p_critical; % exit pressure equals critical pressure if at chocked flow, using equation 18
              temp_exit = (2 / (gamma + 1)) * temp_air_thrust_stage; % exit temperature, using equation 18
              rho_exit = p_exit / (gas_const_air * temp_exit); % exit density, using equation 18

              vel_exit = sqrt(gamma * gas_const_air * temp_exit); % exit velocity, using equation 17

            else % if not chocked flow (which is always)

                mach_exit = sqrt((2 / (gamma - 1)) * ...
                    (((p_air_thrust_stage / p_amb)^((gamma - 1) / gamma)) - 1)); % exit mach number, using equation 19
                p_exit = p_amb; % exit pressure equals ambient pressure, using equation 20
                temp_exit = temp_air_thrust_stage / (1 + ((gamma - 1) / 2) * ...
                    (mach_exit)^(2)); % exit temperature, using equation 20
                rho_exit = p_exit / (gas_const_air * temp_exit); % exit density, using equation 20

                vel_exit = mach_exit * sqrt(gamma * gas_const_air * temp_exit); % exit velocity, using equation 21

            end   
        
            dmassair_dt = -1 * ...
                (discharge_coeff * rho_exit * area_throat * vel_exit); % change in mass of the air in the rocket over change in time, using equation 23

            thrust(i) = (-1 * dmassair_dt * vel_exit) + (p_exit - p_amb) * ...
                area_throat; % thrust of the rocket, using equation 22

        else % if the bottle rocket is in stage 3 flight, no thrust
        
            thrust(i) = 0; % no thrust in stage 3
            
            if phase2_end_i == 0
                
                phase2_end_i = i-1; % finds the indice for the end of phase 1
                
            end
            
        end
        
    end

    % return the thrust values and phase ending indice values
    thrustvec = thrust;
    phase1_end = phase1_end_i;
    phase2_end = phase2_end_i;
    
end