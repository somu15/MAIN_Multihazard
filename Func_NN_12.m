function [req] = Func_NN_12(x,tr,pr,cdf,info,t)

phi22 = interp1(info(1,:),info(2,:),t-x,'linear','extrap');

fd = interp1(tr,pr(1,:),x,'linear','extrap');

pdf = Differentiation(cdf(1,2)-cdf(1,1),cdf(2,:));
phi12 = interp1(cdf(1,:),pdf,x,'linear','extrap');

req = fd.*phi22.*phi12;

end