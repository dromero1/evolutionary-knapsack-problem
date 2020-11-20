function [found,x_star,R_delta,z_delta,X_lnd,Z_lnd] = kp_vnd_1add(x0,n,W,A,R,b,z,ndf)
%KP_VND_1ADD KP Variable neighborhood descent 1-add
%
%   Inputs:
%   x0 - Solution
%   n - Number of items
%   W - Objective coefficients
%   A - Constraint coefficients
%   R - Current resource consumption
%   b - Resource capacity
%   z - Current objective values
%   ndf - Look for local non-dominated solutions
%
%   Outputs:
%   found - A better solution was found
%   x_star - Improved solution
%   R_delta - Change in resource consumption
%   z_delta - Change in objective values
%   X_lnd - Local non-dominated solution
%   Z_lnd - Objetive values of local non-dominated solutions

% Current solution
x = x0;

% Best solution so far
x_star = x;
z_star = z;

% Deltas
R_delta = b*0;
z_delta = z*0;

% Best move so far
b_move = 0;

% Better solution found flag
found = false;

% Local non-dominated solutions
X_lnd = [];
Z_lnd = [];

% Explore neighborhood
for i = 1:n
    if x(i) == false
        % Determine if movement is legal
        if R + A(:,i) <= b
            % New solution objetive values
            z_prime = z + W(:,i);
            % Determine if new solution dominates the
            % current solution
            if prod(z_prime>=z) == 1 && sum(z_prime>z) >= 1
                % Update found flag
                found = true;
                % Determine if new solution dominates the
                % best solution
                if prod(z_prime>=z_star) == 1 && sum(z_prime>z_star) >= 1
                    % Update solution
                    b_move = i;
                    % Update best move objective values
                    z_star = z_prime;
                end
                % Determine if the new solution is not
                % dominated by the current solution
            elseif ~(prod(z>=z_prime) == 1 && sum(z>z_prime) >= 1) && ndf == true
                % Apply movement
                x_nd = x;
                x_nd(i) = ~x_nd(i);
                % Save local non-dominated solution
                X_lnd = [X_lnd; x_nd'];
                Z_lnd = [Z_lnd; z_prime'];
                % Get non-dominated solutions
                [ND,~] = pareto_dominance(Z_lnd);
                X_lnd = X_lnd(ND,:);
                Z_lnd = Z_lnd(ND,:);
            end
        end
    end
end

if found == true
    % Get best movement
    i = b_move;
    % Update best solution
    x(i) = true;
    x_star = x;
    % Update deltas
    R_delta = A(:,i);
    z_delta = W(:,i);
end

end