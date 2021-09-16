function [Xshow] = discrete(X,Prescision)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Xshow = round(X./Prescision).*Prescision;
end

