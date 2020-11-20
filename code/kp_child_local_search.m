function x = kp_child_local_search(x,n,W,A,b,t0,mt)
%KP_CHILD_LOCAL_SEARCH Local search for children individuals
%
%   Inputs:
%   x - Solution
%   n - Number of items
%   W - Objective coefficients
%   A - Constraint coefficients
%   b - Resource capacity
%   t0 - Initial time
%   mt - Maximum execution time
%
%   Outputs:
%   x - Mutated solution

% Current resource consumption
R = A*x';

% Objective values
z = W*x';

% Random neighborhood
j = randi(3);

% Neighborhood-based local search
if j == 1
    [found,x_star,~,~,~,~] = kp_vnd_1add(x,n,W,A,R,b,z,false,t0,mt);
elseif j == 2
    [found,x_star,~,~,~,~] = kp_vnd_2flip(x,n,W,A,R,b,z,false,t0,mt);
else
    [found,x_star,~,~,~,~] = kp_vnd_21flip(x,n,W,A,R,b,z,false,t0,mt);
end

if found == true
    x = x_star;
end

end