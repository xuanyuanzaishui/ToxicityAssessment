%% 绘制用作种群的AP
clear
%% SR
%
load('TestPop\SR\SRpop.mat') % population data
arrhy = 1;
for i = 1:length(SRBaseCells20141)
    ti = SRBaseCells20141(i).times;
    Vi = SRBaseCells20141(i).V;
    Caii = SRBaseCells20141(i).Cai;

    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP(arrhy).times = t-1000;
    cellAP(arrhy).V = V;
    arrhy = arrhy + 1;
end

% 绘图
c1 = [0, 112/255, 192/255];
figure
for i = 1:4:length(cellAP)

    p1 = plot(cellAP(i).times,cellAP(i).V,'linewidth',0.5,'Color',c1);

    set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
    'YLim',[-100 50],'YTick',[-100 -50 0 50],'TickDir','out','XTick',[0 100 200 1000])
hold on
end
% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}



%% 数据拆分
function [times,volts,cais] = splitdata(Vi,Caii,Ti)   
            PCL = 1000;
            i=fix(length(Ti)./PCL);
            intervals = find(~mod(Ti,PCL));%查找ti取模为0的地方
            times = Ti(intervals(i):intervals(i+1),1);%分离已模拟的细胞的时间数据，用以后续检查数据的正确与否
            volts = Vi(intervals(i):intervals(i+1),1);%分离已模拟的细胞的电压数据
            cais = Caii(intervals(i):intervals(i+1),1);
           
end
