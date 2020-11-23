function ch = kp_crossover(p1,p2,rog,n,m,W,A,b,alpha)
%KP_Crossover Crossover operator
%
%   Inputs:
%   p1 - First parent
%   p2 - Second parent
%   rog - Random offspring generation threshold
%   n - Number of items
%   m - Number of constraints
%   W - Objective coefficients
%   A - Constraint coefficients
%   b - Resource capacity
%   alpha - Best candidate percentage
%
%   Outputs:
%   ch - Child

% Child
ch = false(1,n);

% Hamming distance
hdist = pdist([double(p1); double(p2)],'hamming');

if hdist > rog
    % Crossover point
    cp = randi(n-1);
    % Single-point crossover
    ch(1:cp) = p1(1:cp);
    ch(cp+1:end) = p2(cp+1:end);
else
    % Random offspring generation
    ch_raw = kp_grasp_construct_solution(n,m,W,A,b,alpha);
    ch = ch_raw';
end

end