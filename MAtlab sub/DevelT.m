function [A_matrix] =DevelT(GMass,GDiff,A_matrix,bc,not,dt,GVect)
%%development over time
theta =1;
for i=1:not %timesteps
    KDiff_star = GMass/dt +theta*GDiff;    %respect time
    fVec_star = (GMass/dt -GDiff*(1-theta))*GVect;
    an1 = Solve2(KDiff_star,fVec_star,bc)
    A_matrix(:,i+1)  = an1;
    an =an1

end