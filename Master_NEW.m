% %% Single hazard
% 
% clr
% 
% IM = 65;
% 
% ld1 = 0.01; ld2 = 0.01;
% % l12 = 0.01; l23 = 0.05;
% 
% t = 1:10:2000;
% 
% tr = 1:1:2000;
% pr(1,:) = expcdf(tr,1/ld1);
% pr(2,:) = expcdf(tr,1/ld2);
% 
% [tot_time] = Rep_dists_Li_Ell_Hurr_Char_SMPRESS(tr, IM);
% 
% cdf23 = [tr;tot_time(2,:)];
% cdf12 = [tr;tot_time(1,:)];
% 
% for ii = 1:length(t)
%     
%     r23(ii) = quadgk(@(x)Func_NN_23(x,tr,pr,cdf23),0,t(ii));
%     
% end
% % fdd = expcdf(tr,1/l);
% % f = fdd.*exppdf(tr,1./l23);
% % c23 = trapz(tr,f);
% 
% fd = expcdf(t,1/ld1);
% fdd = expcdf(t,1/ld2);
% 
% c23 = quadgk(@(x)Func_CONST(x,t,fdd,cdf23),0,2000);
% r23 = r23/c23;
% 
% 
% 
% 
% 
% Phi23 = interp1(cdf23(1,:),cdf23(2,:),t,'linear','extrap');
% 
% r22 = 1-fdd.*Phi23;
% 
% for ii = 1:length(t)
%     
%     r12(ii) = quadgk(@(x)Func_NN_12(x,t,[fd;fdd],cdf12,[t;r22],t(ii)),0,t(ii));
%     r13(ii) = quadgk(@(x)Func_NN_13(x,t,[fd;fdd],cdf12,[t;r23],t(ii)),0,t(ii));
%     
% end
% 
% c12 = quadgk(@(x)Func_CONST(x,t,fd,cdf12),0,2000);
% 
% r12 = r12/c12;
% r13 = r13/c12;
% 
% 
% 
% Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
% 
% r11 = 1-fd.*Phi12;
% 
% c = r11+r12+r13;
% 
% rec = 0.5*r12+r13;
% plot(t,rec./c)
% 
% 
% % r1s = exp(-l12*t);
% % r2s = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
% % r3s = 1-r1s-r2s;
% % recs = 0.5*r2s+r3s;
% % 
% % plot(t,rec,t,recs)
% % 
% % ylim([0 1])
% 
% clearvars -except IM
% 
% 
% 
% ld1 = 1000; ld2 = 1000;
% % l12 = 0.01; l23 = 0.05;
% 
% t = 1:10:2000;
% 
% tr = 1:1:2000;
% pr(1,:) = expcdf(tr,1/ld1);
% pr(2,:) = expcdf(tr,1/ld2);
% 
% [tot_time] = Rep_dists_Li_Ell_Hurr_Char_SMPRESS(tr, IM);
% 
% cdf23 = [tr;tot_time(2,:)];
% cdf12 = [tr;tot_time(1,:)];
% 
% for ii = 1:length(t)
%     
%     r23(ii) = quadgk(@(x)Func_NN_23(x,tr,pr,cdf23),0,t(ii));
%     
% end
% 
% 
% % fdd = expcdf(tr,1/l);
% % f = fdd.*exppdf(tr,1./l23);
% % c23 = trapz(tr,f);
% 
% fd = expcdf(t,1/ld1);
% fdd = expcdf(t,1/ld2);
% 
% c23 = quadgk(@(x)Func_CONST(x,t,fdd,cdf23),0,2000);
% r23 = r23/c23;
% 
% c12 = quadgk(@(x)Func_CONST(x,t,fd,cdf12),0,2000);
% 
% Phi23 = interp1(cdf23(1,:),cdf23(2,:),t,'linear','extrap');
% 
% r22 = 1-fdd.*Phi23;
% 
% for ii = 1:length(t)
%     
%     r12(ii) = quadgk(@(x)Func_NN_12(x,t,[fd;fdd],cdf12,[t;r22],t(ii)),0,t(ii));
%     r13(ii) = quadgk(@(x)Func_NN_13(x,t,[fd;fdd],cdf12,[t;r23],t(ii)),0,t(ii));
%     
% end
% r12 = r12/c12;
% r13 = r13/c12;
% 
% Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
% 
% r11 = 1-fd.*Phi12;
% 
% c = r11+r12+r13;
% 
% recs = 0.5*r12+r13;
% 
% hold on
% plot(t,recs./c)

% %% With dependence Exponential
% 
% clr
% 
% sims = 1;
% 
% for mm = 1:sims
% 
%    
%     
% NH = 4;%poissrnd(5);
% % dt = exprnd(300,NH-1,1);
% T = 2000;
% dt = [0 150 800 250 T];%[0;dt;T];
% 
% l12 = 1/250; l23 = 1/300;
% 
% 
% 
% t = 1:10:T;
% tr = 1:1:T;
% 
% pr(1,:) = ones(1,length(tr));
% pr(2,:) = ones(1,length(tr));
% 
% cdf23 = [tr;expcdf(tr,1/l23)];
% cdf12 = [tr;expcdf(tr,1/l12)];
% 
% cumtim = dt(1);
% 
% for ss = 1:NH
% 
% 
% 
% for ii = 1:length(t)
%     
%     r23(ii) = quadgk(@(x)Func_NN_23(x,tr,pr,cdf23),0,t(ii));
%     
% end
% 
% % fd = expcdf(t,1/ld1);
% % 
% 
% c23 = quadgk(@(x)Func_CONST(x,tr,pr(2,:),cdf23),0,T);
% r23 = r23/c23;
% 
% 
% 
% fdd = interp1(tr,pr(2,:),t,'linear','extrap');
% Phi23 = interp1(cdf23(1,:),cdf23(2,:),t,'linear','extrap');
% r22 = 1-fdd.*Phi23;
% 
% for ii = 1:length(t)
%     
%     r12(ii) = quadgk(@(x)Func_NN_12(x,tr,pr,cdf12,[t;r22],t(ii)),0,t(ii));
%     r13(ii) = quadgk(@(x)Func_NN_13(x,tr,pr,cdf12,[t;r23],t(ii)),0,t(ii));
%     
% end
% 
% c12 = quadgk(@(x)Func_CONST(x,tr,pr(1,:),cdf12),0,T);
% 
% r12 = r12/c12;
% r13 = r13/c12;
% 
% 
% fd = interp1(tr,pr(1,:),t,'linear','extrap');
% Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
% r11 = 1-fd.*Phi12;
% 
% c = r11+r12+r13;
% 
% rec(ss,:) = (0.5*r12+r13)./c;
% 
% cumtim(ss+1) = cumtim(ss)+dt(ss+1);
% 
% pr(1,:) = interp1(t,1-r11,dt(ss+1)+tr,'linear','extrap');
% pr(2,:) = interp1(t,1-r22,dt(ss+1)+tr,'linear','extrap');
% 
% end
% % plot(t,rec./c)
% 
% % r1s = exp(-l12*t);
% % r2s = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
% % r3s = 1-r1s-r2s;
% % recs = 0.5*r2s+r3s;
% 
% % hold on
% % plot(t,recs)
% 
% % clearvars
% ind1 = 1;
% ind2 = 0;
% for ii = 1:NH
%     
%    tf = t(1):dt(ii+1);
%    ind2 = ind2+length(tf);
%    recfw(ind1:ind2) = interp1(t,rec(ii,:),tf,'linear','extrap');
%    ind1 = ind2+1;
%     
% end
% 
% tff = 1:3000;
% recff(mm,:) = interp1(1:length(recfw),recfw,tff,'linear','extrap');
% 
% progressbar(mm/sims)
% 
% end
% %% Without dependence Exponential
%  
% clr
% 
% sims = 1;
% 
% for mm = 1:sims
% 
%    
%     
% NH = 4;%poissrnd(5);
% % dt = exprnd(300,NH-1,1);
% T = 2000;
% dt = [0 150 800 250 T];%[0;dt;T];
% 
% l12 = 1/250; l23 = 1/300;
% 
% 
% 
% t = 1:10:T;
% tr = 1:1:T;
% 
% pr(1,:) = ones(1,length(tr));
% pr(2,:) = ones(1,length(tr));
% 
% cdf23 = [tr;expcdf(tr,1/l23)];
% cdf12 = [tr;expcdf(tr,1/l12)];
% 
% cumtim = dt(1);
% 
% for ss = 1:NH
% 
% 
% 
% for ii = 1:length(t)
%     
%     r23(ii) = quadgk(@(x)Func_NN_23(x,tr,pr,cdf23),0,t(ii));
%     
% end
% 
% % fd = expcdf(t,1/ld1);
% % 
% 
% c23 = quadgk(@(x)Func_CONST(x,tr,pr(2,:),cdf23),0,T);
% r23 = r23/c23;
% 
% 
% 
% fdd = interp1(tr,pr(2,:),t,'linear','extrap');
% Phi23 = interp1(cdf23(1,:),cdf23(2,:),t,'linear','extrap');
% r22 = 1-fdd.*Phi23;
% 
% for ii = 1:length(t)
%     
%     r12(ii) = quadgk(@(x)Func_NN_12(x,tr,pr,cdf12,[t;r22],t(ii)),0,t(ii));
%     r13(ii) = quadgk(@(x)Func_NN_13(x,tr,pr,cdf12,[t;r23],t(ii)),0,t(ii));
%     
% end
% 
% c12 = quadgk(@(x)Func_CONST(x,tr,pr(1,:),cdf12),0,T);
% 
% r12 = r12/c12;
% r13 = r13/c12;
% 
% 
% fd = interp1(tr,pr(1,:),t,'linear','extrap');
% Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
% r11 = 1-fd.*Phi12;
% 
% c = r11+r12+r13;
% 
% rec(ss,:) = (0.5*r12+r13)./c;
% 
% cumtim(ss+1) = cumtim(ss)+dt(ss+1);
% 
% pr(1,:) = ones(1,length(tr));
% pr(2,:) = ones(1,length(tr));
% 
% end
% % plot(t,rec./c)
% 
% % r1s = exp(-l12*t);
% % r2s = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
% % r3s = 1-r1s-r2s;
% % recs = 0.5*r2s+r3s;
% 
% % hold on
% % plot(t,recs)
% 
% % clearvars
% ind1 = 1;
% ind2 = 0;
% for ii = 1:NH
%     
%    tf = t(1):dt(ii+1);
%    ind2 = ind2+length(tf);
%    recfw(ind1:ind2) = interp1(t,rec(ii,:),tf,'linear','extrap');
%    ind1 = ind2+1;
%     
% end
% 
% tff = 1:3000;
% recff(mm,:) = interp1(1:length(recfw),recfw,tff,'linear','extrap');
% 
% progressbar(mm/sims)
% 
% end

%% With dependence Multihazards

clr

IMe = 0.4; IMh = 80;

sims = 100;

for mm = 1:sims

   
    
NH = 2;
dt = exprnd(365,NH-1,1);
T = 2000;
dt = [0;dt;T];

% l12 = 1/250; l23 = 1/300;



t = 1:10:T;
tr = 1:1:T;

pr(1,:) = ones(1,length(tr));
pr(2,:) = ones(1,length(tr));

[tot_time] = Rep_dists_Li_Ell_Hurr_Char_SMPRESS(tr, IMh);

cdf23 = [tr;tot_time(2,:)];
cdf12 = [tr;tot_time(1,:)];

cumtim = dt(1);

for ss = 1:NH



for ii = 1:length(t)
    
    r23(ii) = quadgk(@(x)Func_NN_23(x,tr,pr,cdf23),0,t(ii));
    
end

% fd = expcdf(t,1/ld1);
% 

c23 = quadgk(@(x)Func_CONST(x,tr,pr(2,:),cdf23),0,T);
r23 = r23/c23;



fdd = interp1(tr,pr(2,:),t,'linear','extrap');
Phi23 = interp1(cdf23(1,:),cdf23(2,:),t,'linear','extrap');
r22 = 1-fdd.*Phi23;

for ii = 1:length(t)
    
    r12(ii) = quadgk(@(x)Func_NN_12(x,tr,pr,cdf12,[t;r22],t(ii)),0,t(ii));
    r13(ii) = quadgk(@(x)Func_NN_13(x,tr,pr,cdf12,[t;r23],t(ii)),0,t(ii));
    
end

c12 = quadgk(@(x)Func_CONST(x,tr,pr(1,:),cdf12),0,T);

r12 = r12/c12;
r13 = r13/c12;


fd = interp1(tr,pr(1,:),t,'linear','extrap');
Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
r11 = 1-fd.*Phi12;

c = r11+r12+r13;

rec(ss,:) = (0.5*r12+r13)./c;

cumtim(ss+1) = cumtim(ss)+dt(ss+1);

pr(1,:) = interp1(t,1-r11,dt(ss+1)+tr,'linear','extrap');
pr(2,:) = interp1(t,1-r22,dt(ss+1)+tr,'linear','extrap');

[tot_time] = Rep_dists_Li_Ell_Eq_Char_SMPRESS(tr, IMe);

cdf23 = [tr;tot_time(2,:)];
cdf12 = [tr;tot_time(1,:)];

end
% plot(t,rec./c)

% r1s = exp(-l12*t);
% r2s = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
% r3s = 1-r1s-r2s;
% recs = 0.5*r2s+r3s;

% hold on
% plot(t,recs)

% clearvars
ind1 = 1;
ind2 = 0;
for ii = 1:NH
    
   tf = t(1):dt(ii+1);
   ind2 = ind2+length(tf);
   recfw(ind1:ind2) = interp1(t,rec(ii,:),tf,'linear','extrap');
   ind1 = ind2+1;
    
end

tff = 1:3000;
recff(mm,:) = interp1(1:length(recfw),recfw,tff,'linear','extrap');

ind = find(recff(mm,:)>1);
recff(mm,ind) = 1;

progressbar(mm/sims)

end
%% Without dependence Multihazards
 
clr

IMe = 0.4; IMh = 80;

sims = 101;

for mm = 1:sims

   
    
NH = 2;
dt = exprnd(365,NH-1,1);
T = 2000;
dt = [0;dt;T];

% l12 = 1/250; l23 = 1/300;



t = 1:10:T;
tr = 1:1:T;

pr(1,:) = ones(1,length(tr));
pr(2,:) = ones(1,length(tr));

[tot_time] = Rep_dists_Li_Ell_Hurr_Char_SMPRESS(tr, IMh);

cdf23 = [tr;tot_time(2,:)];
cdf12 = [tr;tot_time(1,:)];

cumtim = dt(1);

for ss = 1:NH



for ii = 1:length(t)
    
    r23(ii) = quadgk(@(x)Func_NN_23(x,tr,pr,cdf23),0,t(ii));
    
end

% fd = expcdf(t,1/ld1);
% 

c23 = quadgk(@(x)Func_CONST(x,tr,pr(2,:),cdf23),0,T);
r23 = r23/c23;



fdd = interp1(tr,pr(2,:),t,'linear','extrap');
Phi23 = interp1(cdf23(1,:),cdf23(2,:),t,'linear','extrap');
r22 = 1-fdd.*Phi23;

for ii = 1:length(t)
    
    r12(ii) = quadgk(@(x)Func_NN_12(x,tr,pr,cdf12,[t;r22],t(ii)),0,t(ii));
    r13(ii) = quadgk(@(x)Func_NN_13(x,tr,pr,cdf12,[t;r23],t(ii)),0,t(ii));
    
end

c12 = quadgk(@(x)Func_CONST(x,tr,pr(1,:),cdf12),0,T);

r12 = r12/c12;
r13 = r13/c12;


fd = interp1(tr,pr(1,:),t,'linear','extrap');
Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
r11 = 1-fd.*Phi12;

c = r11+r12+r13;

rec(ss,:) = (0.5*r12+r13)./c;

cumtim(ss+1) = cumtim(ss)+dt(ss+1);

pr(1,:) = ones(1,length(tr));
pr(2,:) = ones(1,length(tr));

[tot_time] = Rep_dists_Li_Ell_Eq_Char_SMPRESS(tr, IMe);

cdf23 = [tr;tot_time(2,:)];
cdf12 = [tr;tot_time(1,:)];

end
% plot(t,rec./c)

% r1s = exp(-l12*t);
% r2s = l12/(l23-l12)*(exp(-l12*t)-exp(-l23*t));
% r3s = 1-r1s-r2s;
% recs = 0.5*r2s+r3s;

% hold on
% plot(t,recs)

% clearvars
ind1 = 1;
ind2 = 0;
for ii = 1:NH
    
   tf = t(1):dt(ii+1);
   ind2 = ind2+length(tf);
   recfw(ind1:ind2) = interp1(t,rec(ii,:),tf,'linear','extrap');
   ind1 = ind2+1;
    
end

tff = 1:3000;
recff(mm,:) = interp1(1:length(recfw),recfw,tff,'linear','extrap');

ind = find(recff(mm,:)>1);
recff(mm,ind) = 1;

progressbar(mm/sims)

end
