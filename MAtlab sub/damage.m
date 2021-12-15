function [c] = damage(tburn,t,dt,temperature,hol)
JJ = (t-tburn)

for i = JJ:t
    dam(1,i) = 2*10^98*exp(-12017/(temperature(hol,18+i)-273.15))
    intergrate_trap(1,i) = trapz(dam)
    
  
    
end
c= intergrate_trap.*dt