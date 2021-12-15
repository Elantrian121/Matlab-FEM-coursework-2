%%main
xmin = [0,0.00166667,0.005];
xmax = [0.00166667,0.005,0.01];
ne = [3,6,9];
Net = 18; %total elements
NnT = 19; %total number of nodes
Iter = length(ne)
Tin = 310.15;
Tend = 310.15;
%%parameters
rho = 1200; h_cap = 3300;
rho_b = 1060; h_cap_b = 3770;
T_b = 310.15; G = 0.0375;
f_term = ((G*rho_b*h_cap_b)/rho*h_cap);
D=1;
lambda =3;

%%imesh generation 
%%element matrix
dxel = (xmax(1)-xmin(1))/ne(1);
dxetwo = (xmax(2)-xmin(2))/ne(2);
dxthree = (xmax(3)-xmin(3))/ne(3);
lxl = [xmin(1):dxel:xmax(1)];
lxll = [xmin(2):dxetwo:xmax(2)];
lxlll = [xmin(3):dxthree:xmax(3)];
x_arr = [1:1:NnT];%%building array of coordinates
x_arr(1:4) = lxl;
x_arr(4:10) = lxll;
x_arr(10:19) = lxlll;
%%init matrices and vect
K = zeros(NnT,NnT);%empty stiffness matrix
M = zeros(NnT,NnT);%empty mass matrix
F = zeros(NnT,1);%empty vector
Fin = zeros(NnT,1);%inital vector

%assembling matrices
for i = 1:ne(1)
    dx = (xmax(1)-xmin(1))/ne(1)
    D = 20 /(rho*h_cap);
    Ke = Diff_Matrix(dx,D)
    %add into position in global matrix
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+Ke
end
for i=4:9
    dx = (xmax(2)-xmin(2))/(ne(2))
    D = 40/(rho*h_cap);
    Ke = Diff_Matrix(dx,D)
    %add into position in global matrix
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+Ke
end
for i =10:18
    dx = (xmax(3)-xmin(3))/ne(3)
    D = 25/(rho*h_cap);
    Ke = Diff_Matrix(dx,D)
    %add into position in global matrix
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+Ke
end


%%nbow do the same for the reaction terms
for i=1:ne(1)
    dx = (xmax(1)-xmin(1))/ne(1)
    J = dx*0.5; %%as using equally spaced mesh
    lambda = 1/(rho*h_cap);
    reaction = Reaction_elem(lambda,J);
    %insert to the global matrix in positions
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)-reaction;
end
for i=4:9
    dx = (xmax(2)-xmin(2))/(ne(2))
    J = dx*0.5; %%as using equally spaced mesh
    lambda = ((G*rho_b*h_cap_b)/rho*h_cap)
    reaction = Reaction_elem(lambda,J);
    %insert to the global matrix in positions
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)-reaction;
end
for i=10:18
    dx = (xmax(3)-xmin(3))/ne(3)
    J = dx*0.5; %%as using equally spaced mesh
    lambda = ((G*rho_b*h_cap_b)/rho*h_cap)
    reaction = Reaction_elem(lambda,J);
    %insert to the global matrix in positions
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)-reaction;
end
%%now for the mass matrix
for i = 1:ne(1)
    
    Me = mass_elem(xmin(1),xmax(1),ne(1),rho);
    %insert into global mass matrix
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me;
end
for i =4:9
   
    Me = mass_elem(xmin(2),xmax(2),ne(2),rho);
    %insert into global mass matrix
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me;
end
for i=10:18
    
    Me = mass_elem(xmin(3),xmax(3),ne(3),rho);
    %insert into global mass matrix
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me;
end
%%assem vect
for i=1:ne(1)
    dx = (xmax(1)-xmin(1))/ne(1)
    J = dx*0.5;
    Fe = Source_term(f_term,J);
    F(i:i+1,1) = F(i:i+1,1)+Fe;
end
for i=4:9
    dx = (xmax(2)-xmin(2))/ne(2)
    J = dx*0.5;
    Fe = Source_term(f_term,J);
    F(i:i+1,1) = F(i:i+1,1)+Fe;
end
for i=10:18
    dx = (xmax(3)-xmin(3))/ne(3)
    J = dx*0.5;
    Fe = Source_term(f_term,J);
    F(i:i+1,1) = F(i:i+1,1)+Fe;
end
%% initial conditions
Told = zeros(NnT,1);
Told(:,1) = 310.15;
F(1) = 7.5;
F(NnT) = 7.5;
%%time parameter
dt = 0.5;
t = 50;
%%simplify matrices
Mat1 = (M+0.5*dt*K);
Mat2 = (M-0.5*dt*K);

Tnew = Told;
Fold = F;
Fnew = Fold;

%%time step
for i =1:(t/dt)
    temperature(:,1+i) = Tnew; %store temperature data
    Tnew = Mat1\(Mat2*Told+dt*F);   %crank Nicolson

    Told = Tnew;
   
end

I = find(temperature>317.15); %find values where temp 317.15 is exceeded
tburn = (t-((length(I)/NnT)*dt))%how many of the values there are
time_l = [0:dt:t]
plot(time_l,temperature(1,1:end), 'r')
xlabel("node coordinates")
ylabel("Temperature Kelvin")
title("Enforced condtion to prevent burn")
axis([0 0.01 0 320])
