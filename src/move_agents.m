function [P1, E1] = move_agents(P0, Pdot, E0, dt)
%{  
Let "n" be the number of agents.
--- Inputs ---
P0: 1 by n by 2 array holding positional data for all agents at timestep i
Pdot:  1 by n by 2 array holding positional data for all agents at 
          timestep i+1 BEFORE it has been constrained.
E0:  1 by n array holding energy data for all agents at timestep i
dt: timestep (scalar)
--- Outputs ---
P1:  1 by n by 2 array holding positional data for all agents at 
         timestep i+1 AFTER it has been constrained.
E1: 1 by n array holding energy data for all agents at timestep i+1.
%}

    % Check dimensions explicitly to avoid bugs with single agents
    n = size(P0, 2); 

    % assuming that there are no restrictions on movement
    dvec = Pdot * dt;
    P1 = P0 + dvec;
    E1 = E0; 

    % velocity restriction example --------------------------------------------
    MAXD = 0.1;
    if MAXD > 0
        % Calculate magnitude along the 3rd dimension (coordinates)
        dmag = sqrt(sum(dvec(1,:,:).^2, 3));
        
        for i = 1:n
            if dmag(i) > MAXD
                % Scale the displacement vector to exactly MAXD
                scale_factor = MAXD / dmag(i);
                P1(1,i,:) = P0(1,i,:) + dvec(1,i,:) * scale_factor;
            end
        end
    end

    % energy consumption ------------------------------------------------------
    ECPUD = 20;            % energy consumed per unit distance
    dvec_sat = P1 - P0;    % relative displacement (actual moved distance)
    
    % Recalculate magnitude after constraints
    dmag_sat = sqrt(sum(dvec_sat(1,:,:).^2, 3));  
    
    E1 = E0 - ECPUD * dmag_sat;

    for i = 1:n
        if E1(i) <= 0
            E1(i) = 0;
            % If energy is depleted, agent stays at P0 (did not complete move)
            P1(1,i,:) = P0(1,i,:);
        end
    end
end