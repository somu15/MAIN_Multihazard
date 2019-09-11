% % clr
% % l = 0.01;
% % l2 = 0.05;
% % t = 1:1000;
% % s = 1;
% % p12 = (1-exp(-l*t));
% % p12_n = (1-exp(-l2*t)).*p12;
% % 
% % for ii = 4:length(t)
% %     
% %     s1 = 1:1:t(ii);
% %     
% %     int1 = trapz(s1,exp(-l*s1).*(1-exp(-l2*s1)));
% %     
% % %     mul1(ii) = l*
% %     
% %     res(ii) = int1;
% %     
% % end
% % 
% % plot(t,res/83.26,t,p12_n)
% 
% %% New
% clr
% 
% l = 0.01;
% f = exp(-l*(51:1050));
% 
% l12 = 0.01; l23 = 0.05;
% 
% % t = 1:1000;
% % ind = 1;
% % r12 = (1-f(ind))*trapz(t,exppdf(t,1/l12).*(1-f(ind)-(1-f(ind))*expcdf(t,1/l23)));
% % r23 = trapz(t,exppdf(t,1/l12).*expcdf(t,1/l23));
% for ii = 2:1000
% 
% 
%     t = 1:ii;
%     
% ind = 100;
% 
% % r12(ii) = (1-f(ind))*trapz(t,exppdf(t,1/l12).*(1-f(ind)-(1-f(ind))*expcdf(t,1/l23)));
% % 
% % r23(ii) = (1-f(ind))^2*trapz(t,exppdf(t,1/l12).*expcdf(t,1/l23));
% 
% r11(ii) = 1-f(ind)-(1-f(ind))*expcdf(ii,1/l12);
% 
% r12(ii) = (1-f(ind))*trapz(t,exppdf(t,1/l12).*(1-f(ind)-(1-f(ind))*expcdf(ii-t,1/l23)));
% 
% r13(ii) = (1-f(ind))^2*trapz(t,exppdf(t,1/l12).*(expcdf(ii-t,1/l23)));
% 
% end
% 
% plot(t,r11,t,r12,t,r13)
% 
% % plot(1:1000,r12)
% 
% % for ii = 2:length(t)
% %     
% %    
% %     R12(ii) = exp(-l23*t(ii))*trapz(1:t(ii),(1-f(1:ii)).*exp(-l12*(1:t(ii))));
% %     
% %     R13(ii) = (1-exp(-l23*t(ii)))*trapz(1:t(ii),(1-f(1:ii)).*exp(-l12*(1:t(ii))));
% %     
% %     R11(ii) = exp(-l12*t(ii));
% %     
% %    % R12_O = 
% %     
% % end
%     
% 
% 


%% New 1
% clr
% 
% l = 100;
% % f = exp(-l*(51:1050));
% 
% l12 = 0.05; l23 = 0.01;
% ind = 10;
% 
% for ii = (ind+2):(1000+ind)
% 
%     t = (ind+1):(ii+ind);
%     
% f_1 = exp(-l*(ind));
% f_1inv = 1-f_1;
% 
% f = exp(-l*(t));
% f_inv = 1-f;
% 
% W22 = 1-f-f_inv.*expcdf(ii-t,1/l23); 
% r12(ii-ind-1) = trapz(t,f_1inv*W22.*exppdf(t,1/l12));
% 
% 
% % for jj = 1:length(t)
% %     
% %     t1 = (ind+1):(ii+ind);
% %     if size(t1) == 1
% %         temp13(jj) = 0;
% %     else
% %    temp13(jj) = trapz(t1,exppdf(t1,1/l23));
% %     end
% %     
% % end
% 
% temp13 = 1-exp(-l23*((ii+ind)-t));
% 
% r13(ii-ind-1) = trapz(t,f_1inv*exppdf(t,1/l12).*f_inv.*temp13);
% 
% r11(ii-ind-1) = 1-f_1-f_1inv*expcdf((ii+ind),1/l12);
% 
% % for jj = 1:length(t)
% %     
% %    temp13(jj) = trapz(t,exppdf(t(jj),1/l23)*f_1inv*exppdf(t,1/l12).*f_inv);
% %     
% % end
% % r13(ii-ind-1) = trapz(t,temp13);
% 
% % T13 = exppdf(t,1/l12).*f_inv.*exppdf(ii-t,1/l23);
% 
% 
% 
% end
% 
% r13_a = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
% 
% const = 1;%r11+r12+r13;
% 
% plot(r11./const)
% hold on
% plot(r12./const)
% hold on
% plot(r13./const)
% hold on
% plot(r13_a)

%% New 2
clr

l = 100;
% f = exp(-l*(51:1050));

l12 = 0.05; l23 = 0.01;
ind = 10;

tr = (ind+1):1:(1000+ind);

for ii = 2:length(tr)

    t = 1:tr(ii);
    
f_1 = exp(-l*(ind));
f_1inv = 1-f_1;

f = exp(-l*(t));
f_inv = 1-f;

% W22 = 1-f-f_inv.*expcdf(tr(ii)-t,1/l23); 
% r12(ii) = trapz(t,f_1inv*W22.*exppdf(t,1/l12));

r12(ii) = quadgk(@(t)Func12_EXP(t,tr(ii),l,l12,l23,f_1inv),0,tr(ii));

% for jj = 1:length(t)
%     
%     t1 = (ind+1):(ii+ind);
%     if size(t1) == 1
%         temp13(jj) = 0;
%     else
%    temp13(jj) = trapz(t1,exppdf(t1,1/l23));
%     end
%     
% end

% temp13 = 1-exp(-l23*(tr(ii)-t));
% 
% r13(ii) = trapz(t,f_1inv*exppdf(t,1/l12).*f_inv.*temp13);

r13(ii) = quadgk(@(t)Func23_EXP(t,l,l12,f_1inv,l23,tr(ii)),0,tr(ii));
% 
% r11(ii-ind-1) = 1-f_1-f_1inv*expcdf((ii+ind),1/l12);

% for jj = 1:length(t)
%     
%    temp13(jj) = trapz(t,exppdf(t(jj),1/l23)*f_1inv*exppdf(t,1/l12).*f_inv);
%     
% end
% r13(ii-ind-1) = trapz(t,temp13);

% T13 = exppdf(t,1/l12).*f_inv.*exppdf(ii-t,1/l23);



end

r13_a = l12/(l23-l12)*(exp(-l12*tr)-exp(-l23*tr));

const = 1;%r11+r12+r13;


plot(tr,r12)
hold on
plot(tr,r13)
hold on
plot(r13_a)


%% 05/21/19 (1)

clr

l = 0.01;
l12 = 0.05; l23 = 0.03;

t = 1:0.1:500;

for ii = 1:length(t)
   
    r2(ii) = quadgk(@(x)Func_New_12(x,l,l12,l23,t(ii)),0,t(ii));
    r3(ii) = quadgk(@(x)Func_New_23(x,l,l12,l23,t(ii)),0,t(ii));
    
end
    
% r1 = exp(-l12*t);
% 
% r2s = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
% r3s = 1-r1-r2s;
% 
% c = 1./(r1+r2+r3);
%  
% r1 = r1.*c; r2 = r2.*c; r3 = r3.*c;
% 
% rec = 0.5*r2+r3; rec1 = 0.5*r2s+r3s;
% 
% plot(t,rec,t,rec1)


%% 05/21/19 (2)

clr

l = 0.001;
l12 = 0.005; l23 = 0.01;

t = 1:0.1:500;

for ii = 1:length(t)
   
    r2(ii) = quadgk(@(x)Func_New_12(x,l,l12,l23,t(ii)),0,t(ii));
    r3(ii) = quadgk(@(x)Func_New_23(x,l,l12,l23,t(ii)),0,t(ii));
    
end

r1 = 1-(1-exp(-l*t)).*(1-exp(-l12*t));

% c = r1+r2+r3;
% 
% r1 = r1./c; r2 = r2./c; r3 = r3./c;

rec = 0.5*r2+r3;

r1s = exp(-l12*t);
r2s = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
r3s = 1-r1s-r2s;
recs = 0.5*r2s+r3s;

% plot(t,r1,t,r2,t,r3)

plot(t,rec,t,recs)
