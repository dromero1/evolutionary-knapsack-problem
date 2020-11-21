function ch = kp_crossover(x,y)
%KP_Crossover Crossover operator
%
%   Inputs:
%   x - First parent
%   y - Second parent
%
%   Outputs:
%   ch - Child

% Items count
n = length(x);

% Crossover point
cp = randi(n-1);

% Child
ch = false(1,n);

% Crossover
ch(1:cp) = x(1:cp);
ch(cp+1:end) = y(cp+1:end);

end