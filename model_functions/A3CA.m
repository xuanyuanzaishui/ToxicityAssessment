% 数据装载S
load('TestPop\ClasssIII\male\AFpopICsm1.mat')

Y_CA = [];
ca = 1;
for ii = 1:length(Xin) %for each member of the population 
            Xt = Xin(ii).times;
            Xv = Xin(ii).V;
            Xca = Xin(ii).Cai;
            [times,volts,Cais] = splitdata(Xt,Xv,Xca); % saved the last 10 beats, separate each beat
            
            [peakV,peakdex] = max(volts) ;
            tpeak = t(peakdex) ;

            for i = peakdex : length(times)-1
                if volts(i) < volts(i+1) || V(end) > -50
                    Y_CA(ii) = 0;
                    ca = ca + 1;
                end
             
            end
end

rate = (ca-1)/length(Xin);