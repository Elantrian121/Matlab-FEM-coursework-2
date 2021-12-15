function [R] = Reaction_elem(lambda,J)
%%simplification of the reaction element calculation extract constants so
%%make matrix lambda*J*[2/3,1/3,1/3,2/3]
%typical method of making the matrix
R = zeros(2:2);
R(1,1) = lambda*J*2/3;
R(1,2) = lambda*J*1/3;
R(2,1) = R(1,2);
R(2,2) = R(1,1);
%alternative R = lambda*J*[2/3,1/3;1/3,2/3]would also gen the same matrix
end
