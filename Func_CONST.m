function [req] = Func_CONST(x,t,f,cdf)

fint = interp1(t,f,x,'linear','extrap');

pdf = Differentiation(cdf(1,2)-cdf(1,1),cdf(2,:));
phi23 = interp1(cdf(1,:),pdf,x,'linear','extrap');

req = fint.*phi23;

end