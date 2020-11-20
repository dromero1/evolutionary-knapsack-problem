function X = kp_vnd(x0,n,m,W,A,b,J,t0,mt)
%KP_VND KP Variable neighborhood descent
%
%   Inputs:
%   x0 - Solution
%   n - Number of items
%   m - Number of constraints
%   W - Objective coefficients
%   A - Constraint coefficients
%   b - Resource capacity
%   J - Number of neighborhoods
%   t0 - Initial time
%   mt - Maximum execution time
%
%   Outputs:
%   X - Pareto front of local search

% Solutions in pareto front
X = x0';

% Objetives values in pareto front
Z = X*W';

% Solutions to explore
sol = 1;

% Index
idx = 1;

while sol >= 1 && toc - t0 <= mt
    % Current solution
    x = X(idx,:)';
    % Current resource consumption
    R = A*x;
    % Current objective values
    z = W*x;
    % Local non-dominated solutions
    X_lnd = [];
    Z_lnd = [];
    % Variable neighborhood search
    j = 1;
    while j <= J && toc - t0 <= mt
        if j == 1
            [found,x_star,R_delta,z_delta,X_nd,Z_nd] = kp_vnd_1add(x,n,W,A,R,b,z,true);
        elseif j == 2
            [found,x_star,R_delta,z_delta,X_nd,Z_nd] = kp_vnd_2flip(x,n,W,A,R,b,z,true);
        else
            [found,x_star,R_delta,z_delta,X_nd,Z_nd] = kp_vnd_21flip(x,n,W,A,R,b,z,true);
        end
        if found == true
            % Update current solution
            x = x_star;
            % Update resource consumption
            R = R + R_delta;
            % New objective values
            z = z + z_delta;
            % Update neighborhood counter
            j = 1;
        else
            j = j + 1;
        end
        % Add local non-dominated solutions
        X_lnd = [X_lnd; X_nd];
        Z_lnd = [Z_lnd; Z_nd];
        % Remove duplicates
        [X_lnd,ix,~] = unique(X_lnd,'rows');
        Z_lnd = Z_lnd(ix,:);
    end
    % Save improved solution
    X(idx,:) = x';
    Z(idx,:) = z;
    % Update index
    idx = idx + 1;
    % Update explored solutions
    sol = sol - 1;
    % Number of global solutions
    ns = size(X,1);
    % Add solutions to explore
    lnd = size(X_lnd,1);
    for i = 1:lnd
        % Solution to explore
        x_lnd = X_lnd(i,:);
        z_lnd = Z_lnd(i,:);
        % Evaluate dominance
        nd_flag = true;
        for k = 1:ns
            if prod(Z(k,:)>=z_lnd) == 1 && sum(Z(k,:)>z_lnd) >= 1
                nd_flag = false;
                break;
            end
        end
        if ~ismember(x_lnd,X,'rows') && nd_flag == true
            X = [X; x_lnd];
            Z = [Z; z_lnd];
            sol = sol + 1;
        end
    end
end

% Select solutions
X = X(1:idx-1,:);
Z = Z(1:idx-1,:);

% Get non-dominated solutions
[Ipo,~] = pareto_dominance(Z);
X = X(Ipo,:);

end