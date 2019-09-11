clr

load('C:\Users\lakshd5\Dropbox\Hurricane simulations\Sa_Values.mat');
load('C:\Users\lakshd5\Dropbox\Hurricane simulations\Wind_Values_50.mat');

t_tot = 3000;

sims = 155;

count = 0;

for ii = 1:64
    for jj = 1:500
        if gm(ii,jj)>5
            ind_i(ii,jj) = 1;
        else
            ind_i(ii,jj) = 0;
        end
    end
end

ind = find(sum(ind_i)>0);
gm(:,ind) = [];
gm = gm(:,find(gm(1,:)>0.8));
for ii = 1:sims
     Ts = exprnd(365);
     sto_ts(ii) = Ts;
    for jj = 1:64
        
        IMe = gm(jj,ii);
    IMh = Vgmax_mod(jj,ii);
    
   
    res = STRSAF_MH(IMe, IMh, Ts);
    
    rec_wd(ii,jj,:) = res(1,:);
    rec_wod(ii,jj,:) = res(2,:);
    
    progressbar(jj/64,ii/sims)
    end
end