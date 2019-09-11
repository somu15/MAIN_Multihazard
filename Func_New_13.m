function [req] = Func_New_13(x1,x2,tr,pr,l12,l23,t)

fd = interp1(tr,pr(1,:),x1,'linear','extrap');
fdd = interp1(tr,pr(2,:),x2,'linear','extrap');

phi_12 = exppdf(x1,1/l12);
% Phi_23 = expcdf((t-x),1/l23);
phi_23 = exppdf(x2,1/l23);

% req1 = fd.*phi_12.*fdd.*phi_23;

req = fd.*phi_12.*fdd.*phi_23;


%  for kk = 1:length(x)
%     req2_si(kk) = quadgk(@(x1)second_int(x1,tr,pr,l23),0,t-x(kk));
%  end
% req2 = fdd.*phi_12.*req2_si;
% req = req2;


% req = req2;

% function [req] = second_int(x1,tr,pr,l23)
% 
% fdd = interp1(tr,pr(2,:),x1,'linear','extrap');
% phi_23 = exppdf(x1,1/l23);
% 
% req = fdd.*phi_23;