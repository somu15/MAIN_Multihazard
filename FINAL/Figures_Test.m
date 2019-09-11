% clr
% 
% cd('C:\Users\lakshd5\Desktop\Revs')
% 
% load('WD_MH_Test_D.mat')
% plot(tff,mean(recff),'b')
% hold on
% load('WOD_MH_Test_D.mat')
% plot(tff,mean(recff),['b','--'])
% hold on
% load('WD_MH_Test_D1.mat')
% plot(tff,mean(recff),'r')
% hold on
% load('WOD_MH_Test_D1.mat')
% plot(tff,mean(recff),['r','--'])
% xlim([0 1500])

%% Compute res

clr
load('C:\Users\lakshd5\Dropbox\Revs\New_Results\Test1.mat');

IMe = 1.5; IMhv = [0.5 0.65 0.8 0.95 1.1 1.25 1.4 1.55 1.7]; 

dtv = [50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800];

for ii = 1:length(IMhv)
    
    for jj = 1:length(dtv)
        
        t = (dtv(jj)+1):3000;
        
        rec_wd = reshape(recF_1(ii,jj,:),3000,1);
        rec_wd = rec_wd(tff>dtv(jj)+15);
        
        rec_wod = reshape(recF_2(ii,jj,:),3000,1);
        rec_wod = rec_wod(tff>dtv(jj)+15);
        
        t11 = min(find(rec_wd>0.99));
        res_wd = trapz(1:t11,rec_wd(1:t11));
        
        t11 = min(find(rec_wod>0.99));
        res_wod = trapz(1:t11,rec_wod(1:t11));
        
        per(ii,jj) = abs(res_wd-res_wod)/res_wd*100;
        
    end
    
end


