function [freq] = Conv_INT(x,tr,PHI,R,t)

Rr = interp1(tr,R,t-x,'linear','extrap');
Rr(Rr<0) = 0;

phi = Differentiation(tr(2)-tr(1),PHI);
phir = interp1(tr,phi,x,'linear','extrap');
phir(phir<0) = 0;

freq = phir.*Rr;

end