function A = lloyds_adjacency_matrix(xy, rc)
%{
Let "n" be the number of agents.
--- Inputs ---
xy: n by 2 array holding position data for agents
    e.g., xy(j,2) is the y coordinate for agent-j
rc: radius of communication (scalar) used by all agents
--- Outputs ---
A:  n by n adjacency matrix
    e.g., if A(i,j) = 1, then agent-i receives a signal from agent-j
          if A(i,j) = 0, then agent-i receives NO signal from agent-j
%}

    % 1. Get the number of agents (n)
    n = size(xy, 1);

    % 2. Calculate Pairwise Distances
    % We use implicit expansion (or broadcasting) to compute the distance 
    % between every pair of agents simultaneously.
    
    % xy(:,1) is the column of x-coords. 
    % xy(:,1)' is the row of x-coords.
    % Subtracting them creates an n-by-n matrix of differences.
    dx = xy(:,1) - xy(:,1)';
    dy = xy(:,2) - xy(:,2)';
    
    % Euclidean distance formula: sqrt(dx^2 + dy^2)
    distances = sqrt(dx.^2 + dy.^2);

    % 3. Threshold by Radius of Communication (rc)
    % This creates a logical matrix (true/false), we convert it to double (1/0)
    A = double(distances <= rc);

    % OPTIONAL: Handling Self-Loops
    % Currently, A(i,i) = 1 because the distance from an agent to itself is 0.
    % In many consensus algorithms, agents "communicate" with themselves.
    % If your specific algorithm requires 0s on the diagonal, uncomment the line below:
    % A(1:n+1:end) = 0; 

end