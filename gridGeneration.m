function [xs, N] = gridGeneration(dim, Min, Max, dx)

N = zeros(dim,1);

N = (Max-Min)./dx; % # of points on each axis
N(1) = round(N(1))+1;
N(2) = round(N(2))+1;
% N(3) = 100;
N(3) = round(N(3))+1;
N(4) = round(N(4))+1;

% N = [51; 40; 50];
% N = [101; 78; 100];

xs = zeros(N(1), N(2), N(3), N(4), dim);
for k = 1:N(3)
    tmp = linspace(Min(1), Max(1), N(1))'*ones(1,N(2));
    xs(:,:,k,:,1) = repmat(tmp, [1,1,N(4)]);
    tmp = ones(N(1),1)*linspace(Min(2),Max(2),N(2));
    xs(:,:,k,:,2) = repmat(tmp, [1,1,N(4)]);
    tmp = ones(N(1),N(2))*(Min(3) + (k-1)*dx(3));
    xs(:,:,k,:,3) = repmat(tmp, [1,1,N(4)]);
    xs(:,:,k,:,4) = linspaceNDim(ones(N(1),N(2))*Min(4), ones(N(1),N(2))*Max(4), N(4));
end