%%transient diffusion problem
%% mesh gen
xmin = 0
xmax = 1; %length of domain
ne = 5; %number of elements
dx = xmax/ne; %distance between nodes
x = [xmin:dx:xmax]; %node coordinates
nn = length(x); %number of nodes

%%element matrix
for i= 2:nn
    e((i-1),:) = [i-1 i];
end
%%init matrices and vect
K = zeros(nn,nn);%empty stiffness matrix
M = zeros(nn,nn);%empty mass matrix
F = zeros(nn,1);%empty vector
Fin = zeros(nn,1);%inital vector
%%parameters
D = 1;
lambda = 1;
f_term = 0;
%assembling matrices
for i = 1:ne
    %%kept the difference outside the function as could not seem to pass
    %%the stored variable through into the Diff_matrix function, but
    %%achieves the same.
    Ke = Diff_Matrix(dx,D)
    %add into position in global matrix
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+Ke
end
%%nbow do the same for the reaction terms
for i=1:ne
    J = dx*0.5; %%as using equally spaced mesh
    reaction = Reaction_elem(lambda,J);
    %insert to the global matrix in positions
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+reaction;
end
%%now for the mass matrix
for i = 1:ne
    Me = mass_elem(xmin,xmax,nn);
    %insert into global mass matrix
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me;
end
%%assem vect
for i=1:ne
    J = dx*0.5;
    Fe = Source_term(f_term,J);
    F(i:i+1,1) = F(i:i+1,1)+Fe;
end

%%initial conditions
Told = zeros(nn,1);
F(1) = 0;
F(nn) = 1;
%%time parameters
dt = 0.05;
t = 1;
%%simplify matrices
Mat1 = (M+0.5*dt*K);
Mat2 = (M-0.5*dt*K);
%
Tnew = Told;
Fold = F;
Fnew = Fold;

%%time step
for i =1:(t/dt)
    temperature(:,1+i) = Tnew; %store temperature data
    Tnew = Mat1\(Mat2*Told+dt*F);   %crank Nicolson

    Told = Tnew;
end
%%exact solution
t = 0:dt:1
H = exactSol(dt)
hh = length(H); xx = [1/hh:1/hh:1] 
ex1 = plot(xx,H); Ex = "exact sol"; xlabel("time"), ylabel("solution")
title("exact solution @ x = 0.8")
x = [xmin:dx:xmax]
figure(2)
dtt = length(temperature) 
tili = [1/dtt:1/dtt:1]
a1 = plot(tili,temperature(5,:),'ro'); M1 = "Calculated result";title("x=0.8"); xlabel("time");ylabel("solution");
a2 = plot(x,temperature(:,3),'ro'); M2 = "t = 0.1";title("t=0.1"); xlabel("times");ylabel("solution");
a3 = plot(tili,temperature(5,7),'ro'); M3 = "t = 0.3";title("t=0.3"); xlabel("time");ylabel("solution");
a4 = plot(tili,temperature(5,21),'ro'); M4 = "t = 1.0";title("t=1.0"); xlabel("time");ylabel("solution");
%result comparision
a1 = plot(xx,temperature(5,:),'ro'); M1 = "Calculated result";hold on
ex1 = plot(xx,H); Ex = "Exact sol";title("solution Comparison"); xlabel("time"); ylabel("Temperature")
legend([a1,ex1],[M1,Ex])
hold off
%%L2 norm or root mean square error
for i = 1:hh
    L2(1,i) = rootmean_sq(H(1,i),temperature(5,i))
end
plot(xx,L2,'ro'); title("L2 error @ x= 0.8"); xlabel("time");ylabel("error")

