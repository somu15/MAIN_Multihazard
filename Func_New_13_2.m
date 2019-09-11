function [req] = Func_New_13_2(x,tr,pr,l12,l23,t)

fd = interp1(tr,pr(1,:),x,'linear','extrap');
fdd = interp1(tr,pr(2,:),t-x,'linear','extrap');

phi_12 = exppdf(x,1/l12);
phi_23 = exppdf((t-x),1/l23);

req = fd.*phi_12.*fdd.*phi_23;