function [ft,t]=nilt_STRSAF_paper_REV(F,tm, IM, S1, S2, Ts)
alfa=0; % 0 for exp
M=512; P=2; %2 for exp
N=2*M; wyn=2*P+1;
t=linspace(0,tm,M);
NT=2*tm*N/(N-2); omega=2*pi/NT;
c=alfa+25/NT; s=c-i*omega*(0:N+wyn-2);
Fsc=feval(F, s, IM, S1, S2, tm, Ts);
ft=fft(Fsc(1:N)); ft=ft(1:M);
for n=N:N+wyn-2
 ft(n-N+2,:)=Fsc(n+1)*exp(-i*n*omega*t);
end
ft1=cumsum(ft); ft2=zeros(wyn-1,M);
for I=1:wyn-2
ft=ft2+1./diff(ft1);
 ft2=ft1(2:wyn-I,:); ft1=ft;
end
ft=ft2+1./diff(ft1); ft=2*real(ft)-Fsc(1);
ft=exp(c*t)/NT.*ft; ft(1)=2*ft(1);


function [func_req1] = rec_Hurr1_Ell(s, IM, S1, S2,t_tot, Ts)
        
% t_tot = 2000;

dt = 0.25;
t = dt:dt:t_tot;
% Sa = 2;
% [tot_time] = Rep_dists_EQScaling_SDOF(t, Sa);
% S1 = 1; S2 = 1;
[tot_time] = Rep_dists_Li_Ell_Hurr_Char_SMPRESS(t, IM);
tot_time1 = tot_time(1,:);
tot_time2 = tot_time(2,:);
if length(S1) == 1 && length(S2) == 1

t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = S1*abs(Differentiation(dt,tot_time1));
%tot_time2_diff = Differentiation(dt,tot_time2);
S22 = 0;
else
    ind = find(S1(1,:)<1e-3);
    S1(1,ind) = 0;
    ind = find(S2(1,:)<1e-3);
    S2(1,ind) = 0;
    S11 = interp1(S1(2,:),S1(1,:),t+Ts,'linear','extrap');
    S22 = interp1(S2(2,:),S2(1,:),t+Ts,'linear','extrap');
tot_time1 = (1-S11).*tot_time1;
tot_time2 = (1-S22).*tot_time2;
t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = abs(Differentiation(t(3)-t(2),tot_time1));%[0 (1-S11)].*exppdf(t,1/L1);
S22 = [0 S22];
end

for ii = 1:length(s)
    s1 = s(ii);
        %G12 = LaplTransNUM(t, tot_time1, s);
        G23(ii) = LaplTransNUM(t,[S22]+ tot_time2, s1);
        g12(ii) = LaplTransNUM(t, tot_time1_diff, s1);
        %g23 = LaplTransNUM(t, tot_time2_diff, s);
end
        
        func_req1 = g12.*(1./s-G23);
        
        
function [func_req2] = rec_Hurr2_Ell(s, IM, S1, S2, t_tot, Ts)
        

% t_tot = 2000;

dt = 0.25;
t = dt:dt:t_tot;
% Sa = 2;
% [tot_time] = Rep_dists_EQScaling_SDOF(t, Sa);
% S1 = 1; S2 = 1;

[tot_time] = Rep_dists_Li_Ell_Hurr_Char_SMPRESS(t, IM);
tot_time1 = tot_time(1,:);
tot_time2 = tot_time(2,:);

if length(S1) == 1 && length(S2) == 1
    
t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = S1*abs(Differentiation(dt,tot_time1));
tot_time2_diff = S2*abs(Differentiation(dt,tot_time2));
else
    ind = find(S1(1,:)<1e-3);
    S1(1,ind) = 0;
    ind = find(S2(1,:)<1e-3);
    S2(1,ind) = 0;
    S11 = interp1(S1(2,:),S1(1,:),t+Ts,'linear','extrap');
    S22 = interp1(S2(2,:),S2(1,:),t+Ts,'linear','extrap');
    tot_time1 = (1-S11).*tot_time1;
tot_time2 = (1-S22).*tot_time2;
t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = abs(Differentiation(t(3)-t(2),tot_time1));%[0 (1-S11)].*exppdf(t,1/L1);
tot_time2_diff = abs(Differentiation(t(3)-t(2),tot_time2));%[0 (1-S22)].*exppdf(t,1/L2);
end

for ii = 1:length(s)
    s1 = s(ii);

        %G12 = LaplTransNUM(t, tot_time1, s);
        %G23 = LaplTransNUM(t, tot_time2, s);
        g12(ii) = LaplTransNUM(t, tot_time1_diff, s1);
        g23(ii) = LaplTransNUM(t, tot_time2_diff, s1);
end
        %func_req1 = g12*(1/s-G23);
        func_req2 = g12.*(g23./s);     
        
        
function [func_req1] = rec_Eq1_Ell(s, IM, S1, S2,t_tot, Ts)
        
% t_tot = 2000;

dt = 0.25;
t = dt:dt:t_tot;
% Sa = 2;
% [tot_time] = Rep_dists_EQScaling_SDOF(t, Sa);
% S1 = 1; S2 = 1;
[tot_time] = Rep_dists_Li_Ell_Eq_Char_SMPRESS(t, IM);
tot_time1 = tot_time(1,:);
tot_time2 = tot_time(2,:);
if length(S1) == 1 && length(S2) == 1

t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = S1*abs(Differentiation(dt,tot_time1));
%tot_time2_diff = Differentiation(dt,tot_time2);
S22 = 0;
else
    ind = find(S1(1,:)<1e-3);
    S1(1,ind) = 0;
    ind = find(S2(1,:)<1e-3);
    S2(1,ind) = 0;
    S11 = interp1(S1(2,:),S1(1,:),t+Ts,'linear','extrap');
    S22 = interp1(S2(2,:),S2(1,:),t+Ts,'linear','extrap');
tot_time1 = (1-S11).*tot_time1;
tot_time2 = (1-S22).*tot_time2;
t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = abs(Differentiation(t(3)-t(2),tot_time1));%[0 (1-S11)].*exppdf(t,1/L1);
S22 = [0 S22];
end

for ii = 1:length(s)
    s1 = s(ii);
        %G12 = LaplTransNUM(t, tot_time1, s);
        G23(ii) = LaplTransNUM(t,[S22]+ tot_time2, s1);
        g12(ii) = LaplTransNUM(t, tot_time1_diff, s1);
        %g23 = LaplTransNUM(t, tot_time2_diff, s);
end
        
        func_req1 = g12.*(1./s-G23);
        
        
function [func_req2] = rec_Eq2_Ell(s, IM, S1, S2, t_tot, Ts)
        

% t_tot = 2000;

dt = 0.25;
t = dt:dt:t_tot;
% Sa = 2;
% [tot_time] = Rep_dists_EQScaling_SDOF(t, Sa);
% S1 = 1; S2 = 1;

[tot_time] = Rep_dists_Li_Ell_Eq_Char_SMPRESS(t, IM);
tot_time1 = tot_time(1,:);
tot_time2 = tot_time(2,:);

if length(S1) == 1 && length(S2) == 1
    
t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = S1*abs(Differentiation(dt,tot_time1));
tot_time2_diff = S2*abs(Differentiation(dt,tot_time2));
else
    ind = find(S1(1,:)<1e-3);
    S1(1,ind) = 0;
    ind = find(S2(1,:)<1e-3);
    S2(1,ind) = 0;
    S11 = interp1(S1(2,:),S1(1,:),t+Ts,'linear','extrap');
    S22 = interp1(S2(2,:),S2(1,:),t+Ts,'linear','extrap');
    tot_time1 = (1-S11).*tot_time1;
tot_time2 = (1-S22).*tot_time2;
t = [0 t];
tot_time1 = [0 tot_time1];
tot_time2 = [0 tot_time2];
tot_time1_diff = abs(Differentiation(t(3)-t(2),tot_time1));%[0 (1-S11)].*exppdf(t,1/L1);
tot_time2_diff = abs(Differentiation(t(3)-t(2),tot_time2));%[0 (1-S22)].*exppdf(t,1/L2);
end

for ii = 1:length(s)
    s1 = s(ii);

        %G12 = LaplTransNUM(t, tot_time1, s);
        %G23 = LaplTransNUM(t, tot_time2, s);
        g12(ii) = LaplTransNUM(t, tot_time1_diff, s1);
        g23(ii) = LaplTransNUM(t, tot_time2_diff, s1);
end
        %func_req1 = g12*(1/s-G23);
        func_req2 = g12.*(g23./s);           