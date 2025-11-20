function obset_data = observation_set(p, x, y, d, ro)
%{ 
Let "n" be the number of agents in the observation set.
Let "m' be the length of the arena
--- Inputs ---
p: n by 2 array holding position data for agents in obs set
x, y, d: m by m meshgrids holding arena data 
ro: radius of observation (scalar) used by all agents
--- Outputs ---
obset_data: w by 3 array.
    Column 1: x coordinate weighted by density (x * d)
    Column 2: y coordinate weighted by density (y * d)
    Column 3: density (d)
    
    *IMPORTANT*: We include Column 3 (d) so you can calculate the 
    total mass (sum(d)) required for centroid normalization.
%}

    n = size(p, 1);    % number of agents in observation set
    
    % 1. Create a Boolean Mask for the whole grid
    %    Start with all false (no points observed)
    observed_mask = false(size(x));
    
    % 2. Update mask for each agent
    for i = 1:n
        % Calculate distance from agent i to ALL grid points
        dist_sq = (x - p(i,1)).^2 + (y - p(i,2)).^2;
        
        % Update the mask: points already seen OR points seen by this agent
        observed_mask = observed_mask | (dist_sq <= ro^2);
    end
    
    % 3. Extract data using the combined mask
    %    This automatically handles duplicates/overlaps and is much faster
    %    than unique().
    
    x_vals = x(observed_mask);
    y_vals = y(observed_mask);
    d_vals = d(observed_mask);
    
    % 4. Construct Output
    %    Returns [x*rho, y*rho, rho]
    obset_data = [x_vals .* d_vals, y_vals .* d_vals, d_vals];

end