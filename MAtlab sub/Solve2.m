function [Q] = Solve2(k,f,bc)
[ngn,ngn] = size(k);
fdof = [1:ngn];
%
d = zeros(size(fdof));
Q = zeros(size(fdof));
%boundary condition storing
pdof = bc(:,1); %position
fdof(pdof) = [];
dp = zeros(size(fdof)); dp = dp';
dp(1,1) = bc(1,2); 
dp(length(fdof),1) = bc(2,2)
dpp = bc(:,2); %values
%dp = bc(:,2); %values

G = k(fdof,fdof)

s = G\(f(fdof) - G*dp)
    
d(pdof) = dpp;
d(fdof) = s;
Q = k*d'-f
end
