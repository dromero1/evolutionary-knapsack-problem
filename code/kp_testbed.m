function [rT,rN,rS,rF,rFS,rI1,rI2,rDU] = kp_testbed(s)
%KP_TESTBED Knapsack problem testbed
%
%   Inputs:
%   s - Scenario
%
%   Outputs:
%   rT - Processing times
%   rN - Number of solutions
%   rS - Number of non-dominated solutions
%   rF - Number of non-dominated feasible solutions
%   rFS - Feasible share
%   rI1 - Probability of generating solutions in the Pareto Front
%   rI2 - Probability of generating non-dominated solutions
%   rDU - Mean perc. distance to upper bound

rng('default');

%% Technical parameters

% Instance count
IC = 20;

% Filenames
input_file = 'files/Input.xlsx';
output_file = 'files/Output.xlsx';

% Debug mode
dbg = false;

% Maximum execution time
mt = 10;

%% Results
rT = zeros(IC,1);
rN = zeros(IC,1);
rS = zeros(IC,1);
rF = zeros(IC,1);
rFS = zeros(IC,1);
rI1 = zeros(IC,1);
rI2 = zeros(IC,1);
rDU = zeros(IC,1);

%% Main loop
for ti = 1:IC
    %% Extraction
    % Raw problem
    P_raw = readmatrix(input_file,'Sheet',['I',num2str(ti)],'NumHeaderLines',0);
    % Control record
    cr = P_raw(1,1:3);
    % Items count
    n = cr(1);
    % Restrictions count
    m = cr(2);
    % Objectives count
    p = cr(3);
    % Constraints coefficients
    A = P_raw(2:m+1,1:n);
    % Resources capacity
    b = P_raw(2:m+1,n+1);
    % Objective coefficients
    W = P_raw(m+2:m+p+1,1:n);
    %% Execution
    if s == 1
        MR = kp_scenario1(ti,n,p,m,W,A,b,mt,dbg);
    end
    %% Pareto front
    % Mix solutions
    X = [];
    Z = [];
    for j = 1:length(MR)
        Z = [Z; MR(j,:).Z];
        X = [X; MR(j,:).X];
    end
    % Remove duplicates
    [X,ix,~] = unique(X,'rows');
    Z = Z(ix,:);
    % Remove infeasible solutions
    If = (Z(:,p+1)==1);
    fc = sum(If);
    if fc >= 1
        X = X(If,:);
        Z = Z(If,:);
    end
    F = Z(:,p+1);
    Z = Z(:,1:p);
    % Get non-dominated solutions
    [Ipo,~] = pareto_dominance(Z);
    PX = X(Ipo,:);
    PF = F(Ipo);
    cP = size(PX,1);
    % Upper bound
    ub = abs(W)*ones(n,1);
    % Instance statistics
    fprintf('Instance %d - Statistics\n',ti);
    fprintf('Number of solutions: %d\n',size(X,1));
    fprintf('Number of feasible solutions: %d\n',fc);
    fprintf('Number of pareto-optimal solutions: %d\n',size(PX,1));
    fprintf('Number of feasible solutions in pareto front: %d\n',sum(PF==1));
    for j = 1:length(MR)
        % Method's results
        mr = MR(j,:);
        % Solution count
        cA = size(mr.Z,1);
        cF = sum(mr.Z(:,p+1)==1);
        % Solutions in Pareto Front
        IAP = intersect(mr.X,PX,'rows');
        cAP = size(IAP,1);
        % Calculate statistics
        i1 = cAP / cP;
        i2 = cAP / cA;
        % Distance to upper bound
        Z_prime = mr.Z';
        Z_prime = Z_prime(1:p,:);
        d2ub = (ub-Z_prime)./ub;
        md2ub = mean(d2ub,'all');
        % Save statistics
        rT(ti,mr.mid) = mr.t;
        rN(ti,mr.mid) = mr.nsol;
        rS(ti,mr.mid) = cA;
        rF(ti,mr.mid) = cF;
        rFS(ti,mr.mid) = cF / cA;
        rI1(ti,mr.mid) = i1;
        rI2(ti,mr.mid) = i2;
        rDU(ti,mr.mid) = md2ub;
        % Display
        fprintf('Method %s (time = %0.4f, sol. %d, psol. = %d, ',mr.mtd,mr.t,mr.nsol,cA);
        fprintf('pfea. = %d, md2ub = %0.2f, ',cF,md2ub)
        fprintf('I1 = %0.2f, I2 = %0.2f)\n',i1,i2);
    end
    %% Download results
    nsol = size(PX,1);
    % Build output matrix
    O = NaN(nsol+1,1+max(sum(PX,2))+m+p);
    O(1,1) = nsol;
    for j = 1:nsol
        x = PX(j,:);
        sline = [sum(x) find(x) (A*x')' (W*x')'];
        O(j+1,1:length(sline)) = sline;
    end
    % Write solutions to output file
    writematrix(O,output_file,'Sheet',['I',num2str(ti)]);
end

end