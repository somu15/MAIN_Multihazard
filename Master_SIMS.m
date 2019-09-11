clr

% IMe = 0.4; IMh = 80;

IMe = 1.5; IMhv = 60;%[50 65 80 95 110 125 140 155];

dtv = 150;%[50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800];

sims = 1;
count = 0;
for imi = 1:length(IMhv)
    
    IMh = IMhv(imi);
    for dtt = 1:length(dtv)

        dt = dtv(dtt);
for mm = 1:sims

   
    
NH = 2;
% dt = 150;%exprnd(365,NH-1,1);
T = 2000;
dt = [0;dt;T];

% l12 = 1/250; l23 = 1/300;



t = 1:10:T;
tr = 1:1:T;

pr(1,:) = ones(1,length(tr));
pr(2,:) = ones(1,length(tr));

[tot_time] = Rep_dists_Li_Ell_Eq_Char_SMPRESS(tr, IMe);

cdf23 = [tr;tot_time(2,:)];
cdf12 = [tr;tot_time(1,:)];

cumtim = dt(1);

for ss = 1:NH

if ss == 1

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
    r11(ii) = quadgk(@(x)Func_CONST(x,tr,pr(1,:),cdf12),0,t(ii));
end

c12 = quadgk(@(x)Func_CONST(x,tr,pr(1,:),cdf12),0,T);

r12 = r12/c12;
r13 = r13/c12;
r11 = 1-1/c12*r11;

fd = interp1(tr,pr(1,:),t,'linear','extrap');
Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
% r11 = 1-fd.*Phi12;

c = 1;%r11+r12+r13;

rec_1(ss,:) = (0.5*r12+r13)./c;
rec_2(ss,:) = (0.5*r12+r13)./c;

cumtim(ss+1) = cumtim(ss)+dt(ss+1);



else

    clearvars pr
    
pr1(1,:) = interp1(t,1-r11,dt(ss)+tr,'linear','extrap');
pr1(2,:) = interp1(t,1-r22,dt(ss)+tr,'linear','extrap');

pr2(1,:) = ones(1,length(tr));
pr2(2,:) = ones(1,length(tr));

[tot_time] = Rep_dists_Li_Ell_Hurr_Char_SMPRESS(tr, IMh);

cdf23 = [tr;tot_time(2,:)];
cdf12 = [tr;tot_time(1,:)];


for ii = 1:length(t)
    
    r23_1(ii) = quadgk(@(x)Func_NN_23(x,tr,pr1,cdf23),0,t(ii));
    r23_2(ii) = quadgk(@(x)Func_NN_23(x,tr,pr2,cdf23),0,t(ii));
    
end

% fd = expcdf(t,1/ld1);
% 

c23_1 = quadgk(@(x)Func_CONST(x,tr,pr1(2,:),cdf23),0,T);
r23_1 = r23_1/c23_1;



fdd_1 = interp1(tr,pr1(2,:),t,'linear','extrap');
fdd_2 = interp1(tr,pr2(2,:),t,'linear','extrap');
Phi23 = interp1(cdf23(1,:),cdf23(2,:),t,'linear','extrap');
r22_1 = 1-r23_1;%fdd_1.*Phi23;
r22_2 = 1-r23_2;%fdd_2.*Phi23;

for ii = 1:length(t)
    
    r12_1(ii) = quadgk(@(x)Func_NN_12(x,tr,pr1,cdf12,[t;r22],t(ii)),0,t(ii));
    r13_1(ii) = quadgk(@(x)Func_NN_13(x,tr,pr1,cdf12,[t;r23],t(ii)),0,t(ii));
    
    r12_2(ii) = quadgk(@(x)Func_NN_12(x,tr,pr2,cdf12,[t;r22],t(ii)),0,t(ii));
    r13_2(ii) = quadgk(@(x)Func_NN_13(x,tr,pr2,cdf12,[t;r23],t(ii)),0,t(ii));
    
end

c12_1 = quadgk(@(x)Func_CONST(x,tr,pr1(1,:),cdf12),0,T);
c12_2 = quadgk(@(x)Func_CONST(x,tr,pr2(1,:),cdf12),0,T);

r12_1 = r12_1/c12_1;
r13_1 = r13_1/c12_1;

r12_2 = r12_2/c12_2;
r13_2 = r13_2/c12_2;


fd_1 = interp1(tr,pr1(1,:),t,'linear','extrap');
Phi12 = interp1(cdf12(1,:),cdf12(2,:),t,'linear','extrap');
r11_1 = 1-r12_1-r13_1;%fd_1.*Phi12;

fd_2 = interp1(tr,pr2(1,:),t,'linear','extrap');
r11_2 = 1-r12_2-r13_2;%fd_2.*Phi12;

c_1 = r11_1+r12_1+r13_1; 
c_2 = r11_2+r12_2+r13_2;

rec_1(ss,:) = (0.5*r12_1+r13_1)./c_1;
rec_2(ss,:) = (0.5*r12_2+r13_2)./c_2;

cumtim(ss+1) = cumtim(ss)+dt(ss+1);

end
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
   recfw_1(ind1:ind2) = interp1(t,rec_1(ii,:),tf,'linear','extrap');
   recfw_2(ind1:ind2) = interp1(t,rec_2(ii,:),tf,'linear','extrap');
   ind1 = ind2+1;
    
end

tff = 1:3000;
recff_1 = interp1(1:length(recfw_1),recfw_1,tff,'linear','extrap');
recff_2 = interp1(1:length(recfw_2),recfw_2,tff,'linear','extrap');

ind1 = find(recff_1(mm,:)>1);
ind2 = find(recff_2(mm,:)>1);
recff_1(ind1) = 1;
recff_2(ind2) = 1;



end

recF_1(imi,dtt,:) = recff_1;
recF_2(imi,dtt,:) = recff_2;
count = count + 1;
progressbar(count/(length(IMhv)*length(dtv)))
    end
end

plot(recff_1)
hold on
plot(recff_2)
