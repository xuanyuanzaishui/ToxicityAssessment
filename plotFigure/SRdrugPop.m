clear

load("BasePopCellsEAD.mat",'BaseCellsEAD')
EAD = (1:36);

c2 = [0 119 206]/255;%EAD\fature

%% Amiodarone
%
n = 21;
random_num = EAD(randperm(numel(EAD),n));
random_num = sort(random_num);

% 绘图
figure
for i = 1:length(random_num)
    p2 = plot(BaseCellsEAD(random_num(i)).times,BaseCellsEAD(random_num(i)).V,'linewidth',0.8,'Color',c2);
        
    set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000],'color','none')
hold on
end

% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% Dofetilide
%
n = 22;
random_num = EAD(randperm(numel(EAD),n));
random_num = sort(random_num);

% 绘图
figure
for i = 1:length(random_num)
    p2 = plot(BaseCellsEAD(random_num(i)).times,BaseCellsEAD(random_num(i)).V,'linewidth',0.8,'Color',c2);
        
    set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end

% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% Dronedarone
%
n = 26;
random_num = EAD(randperm(numel(EAD),n));
random_num = sort(random_num);

% 绘图
figure
for i = 1:length(random_num)
    p2 = plot(BaseCellsEAD(random_num(i)).times,BaseCellsEAD(random_num(i)).V,'linewidth',0.8,'Color',c2);
        
    set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end

% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% Ibutilide
%
n = 23;
random_num = EAD(randperm(numel(EAD),n));
random_num = sort(random_num);

% 绘图
figure
for i = 1:length(random_num)
    p2 = plot(BaseCellsEAD(random_num(i)).times,BaseCellsEAD(random_num(i)).V,'linewidth',0.8,'Color',c2);
        
    set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end

% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% Sotalol
%
n = 24;
random_num = EAD(randperm(numel(EAD),n));
random_num = sort(random_num);

% 绘图
figure
for i = 1:length(random_num)
    p2 = plot(BaseCellsEAD(random_num(i)).times,BaseCellsEAD(random_num(i)).V,'linewidth',0.8,'Color',c2);
        
    set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end

% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% Vernakalant
%
n = 25;
random_num = EAD(randperm(numel(EAD),n));
random_num = sort(random_num);

% 绘图
figure
for i = 1:length(random_num)
    p2 = plot(BaseCellsEAD(random_num(i)).times,BaseCellsEAD(random_num(i)).V,'linewidth',0.8,'Color',c2);
        
    set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end

% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% 函数
%数据拆分
function [times,volts,cais] = splitdata(Vi,Caii,Ti)   
            PCL = 1000;
            i=fix(length(Ti)./PCL);
            intervals = find(~mod(Ti,PCL));%查找ti取模为0的地方
            times = Ti(intervals(i):intervals(i+1),1);%分离已模拟的细胞的时间数据，用以后续检查数据的正确与否
            volts = Vi(intervals(i):intervals(i+1),1);%分离已模拟的细胞的电压数据
            cais = Caii(intervals(i):intervals(i+1),1);
           
end
