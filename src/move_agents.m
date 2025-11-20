function [P1, E1] = move_agents(P0, Pdot, E0, dt)
%{  
Let "n" be the number of agents.

--- Inputs ---
P0: 1 by n by 2 array holding positional data for all agents at timestep i
    e.g, P0(1,j,2) is the y position data for agent-j at timestep i 

Pdot:  1 by n by 2 array holding positional data for all agents at 
          timestep i+1 BEFORE it has been constrained. Same structure as
          P0.

E0:  1 by n array holding energy data for all agents at timestep i
    e.g., E0(j) is the energy of agent-j at timestep i 

dt: timestep (scalar)

--- Outputs ---
P1:  1 by n by 2 array holding positional data for all agents at 
         timestep i+1 AFTER it has been constrained. Same structure as
          P0.

E1: 1 by n array holding energy data for all agents at timestep i+1.
    Same structure as E0.
%}

% assuming that there are no restrictions on movement
dvec = Pdot * dt;
P1 = P0 + dvec;

E1 = E0;  % energy can just pass through if you don't want to use it

% velocity restriction example --------------------------------------------
MAXD = 0.1;
if MAXD > 0
    n = length(P0);     % number of agents
    dmag = sqrt(sum(dvec(1,:,:).^2,3));
    for i=1:n
        r = dmag(i) / MAXD;
        if r > 1
            P1(1,i,:) = P0(1,i,:) + dvec(1,i,:) / r;
        end
    end
end

% energy consumption ------------------------------------------------------
ECPUD = 20;            % energy consumed per unit distance
dvec_sat = P1 - P0;     % relative displacement

dmag_sat = sqrt(sum(dvec_sat(1,:,:).^2,3));  
% distance travelled after velocity constrained

E1 = E0 - ECPUD * dmag_sat;
for i = 1:n
    if E1(i) <= 0
        E1(i) = 0;
        P1(1,i,:) = P0(1,i,:);
    end
end

end