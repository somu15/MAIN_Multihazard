function [req] = Func_New_12(x,tr,pr,l12,l23,t)

% f = exp(-l*x);
% finv = 1-f;
% 
% f1 = exp(-l*(t-x));
% f1inv = 1-f1;

% req = l12*exp(-l12*x).*(1-f-finv.*(1-exp(-l23*(t-x))));
% req = finv.*exppdf(x,1/l12).*(1-f1-f1inv.*(1-exp(-l23*(t-x))));

fd = interp1(tr,pr(1,:),x,'linear','extrap');
fdd = interp1(tr,pr(2,:),t-x,'linear','extrap');

phi_12 = exppdf(x,1/l12);
Phi_23 = expcdf((t-x),1/l23);

req = fd.*phi_12.*(1-fdd.*Phi_23);

end