function [P] = evaluateP(Q,beta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Q_max = max(Q);
a = exp(beta*(Q-Q_max));
b = sum(a);

P = a/b;
end

