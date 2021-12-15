function [Me] = mass_elem(xmin,xmax,ne)
nn=ne+1;
mesh.x = linspace(xmin,xmax,nn)
mesh.connect = [1:nn-1;2:nn]
shape = @(x) 0.5*[1-x,1+x]
for conn = mesh.connect
    Me = zeros(2,2)
    le = length(mesh.x(conn))
    for xi = [1.0, -1.0] /sqrt(3)
        N = shape(xi)
        Me = Me + N'*N*le/2
    end
    Me
end