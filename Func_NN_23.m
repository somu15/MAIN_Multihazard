function [req] = Func_NN_23(x,tr,pr,cdf)

fdd = interp1(tr,pr(2,:),x,'linear','extrap');

pdf = Differentiation(cdf(1,2)-cdf(1,1),cdf(2,:));
phi23 = interp1(cdf(1,:),pdf,x,'linear','extrap');

req = fdd.*phi23;


end