function [TotalMesh] = Assemble(p,q,r)
%%generating vales for Supermesh.nvec
TotalMesh.nvec = zeros(1,19);
TotalMesh.nvec(1,1:4) = p;
TotalMesh.nvec(1,4:10) =q;
TotalMesh.nvec(1,10:19) =r;
end
