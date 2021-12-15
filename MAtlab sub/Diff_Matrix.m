function [k] = Diff_Matrix(difference,D)
%%element matrix for diffusion
%below can be achieved with same result in one line
%k=D*(2/(Mesh.elem(i).x(2)-mesh.elem(i).x(1))) 
k = zeros(2:2);
k(1,1) = D*(2/difference);
k(1,2) = -D*(2/difference);
k(2,2) = k(1,1);
k(2,1) = k(1,2);
% now gives diffusion matrix for each element
end
