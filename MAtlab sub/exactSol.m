function [H] = exactSol(dt)
%%loops for time and calls exact solution builds and propagates a matrix
Tend = 1
tStart = 0

t = tStart:dt:Tend
X = zeros(1,length(t))
dx = 1/length(t)
X(:,1:end) = 0:dx:X
H = zeros(1,length(t))
for i=1:length(t)
    H(:,i) = TransientAnalyticSoln(0.8,t(i))
end
