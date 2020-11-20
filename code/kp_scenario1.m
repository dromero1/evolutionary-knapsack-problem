function MR = kp_scenario1(ti,n,p,m,W,A,b,mt,dbg)
%KP_SCENARIO1 Comparing all methods
%
%   Inputs:
%   ti - Test instance
%   n - Number of items
%   p - Number of objectives
%   m - Number of constraints
%   W - Objective coefficients
%   A - Constraint coefficients
%   b - Resource capacity
%   mt - Maximum execution time
%   dbg - Debug mode
%
%   Outputs:
%   MR - Results collection

% Results collection
MR = [];

% Number of neighborhoods
J = 3;

% Method id
mid = 1;

% Alpha sweep
alpha_sweep = [0.05 0.15 0.25];

%% GRASP method
% for alpha = alpha_sweep
%     % Get solutions
%     tic
%     [X,Z,nsol] = kp_grasp(ti,n,p,m,W,A,b,alpha,J,mt,dbg,false);
%     time = toc;
%     % Save results
%     mr.mid = mid;
%     mr.mtd = sprintf('G-%0.2f',alpha);
%     mr.X = X;
%     mr.Z = Z;
%     mr.t = time;
%     mr.nsol = nsol;
%     MR = [MR; mr];
%     % Update method instance id
%     mid = mid + 1;
% end

%% GRASP VND method
% for alpha = alpha_sweep
%     % Get solutions
%     tic
%     [X,Z,nsol] = kp_grasp(ti,n,p,m,W,A,b,alpha,J,mt,dbg,true);
%     time = toc;
%     % Save results
%     mr.mid = mid;
%     mr.mtd = sprintf('G-VND-%0.2f',alpha);
%     mr.X = X;
%     mr.Z = Z;
%     mr.t = time;
%     mr.nsol = nsol;
%     MR = [MR; mr];
%     % Update method instance id
%     mid = mid + 1;
% end

%% GA method
for alpha = alpha_sweep
    % Get solutions
    tic
    [X,Z,nsol] = kp_ga(ti,n,p,m,W,A,b,alpha,J,mt,dbg);
    time = toc;
    % Save results
    mr.mid = mid;
    mr.mtd = sprintf('G-GA-%0.2f',alpha);
    mr.X = X;
    mr.Z = Z;
    mr.t = time;
    mr.nsol = nsol;
    MR = [MR; mr];
    % Update method instance id
    mid = mid + 1;
end

end