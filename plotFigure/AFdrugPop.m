clear

load("BasePopCellsEAD.mat",'BaseCellsEAD')

%%2014 2022
c1 = [184 207 139]/255;
%c2 = [37 109 183]/255;%EAD\fature
c2 = [239 138 67]/255;%EAD\fature
settings.Folder = 'TestPop/AFClassIII/male';%SR
%% Amiodarone
%
settings.SubFolder = 'Vernakalant';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File,'CABaseCells');
load(File1,'Y_Arr')

arrhy = 1;
for i = 1:length(CABaseCells)
    ti = CABaseCells(i).times;
    Vi = CABaseCells(i).V;
    Caii = CABaseCells(i).Cai;
    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP(arrhy).times = t-1000;
    cellAP(arrhy).V = V;
    arrhy = arrhy + 1;
end
% 绘图
figure
for i = 1:length(cellAP)
    if Y_Arr(i) == 0 
        p1 = plot(cellAP(i).times,cellAP(i).V,'linewidth',0.3,'Color',c1);
        
    end
set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end
hold on

p2 = plot(BaseCellsEAD(8).times,BaseCellsEAD(8).V,'linewidth',0.3,'Color',c2);
p12.Color(4) = 0.2;
hold on
p2 = plot(BaseCellsEAD(23).times,BaseCellsEAD(23).V,'linewidth',0.3,'Color',c2);
p12.Color(4) = 0.2;


% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% Dofetilide
%
settings.SubFolder = 'Dofetilide';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File,'CABaseCells');
load(File1,'Y_Arr')

arrhy = 1;
for i = 1:length(CABaseCells)
    ti = CABaseCells(i).times;
    Vi = CABaseCells(i).V;
    Caii = CABaseCells(i).Cai;
    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP(arrhy).times = t-1000;
    cellAP(arrhy).V = V;
    arrhy = arrhy + 1;
end
% 绘图
figure
for i = 1:length(cellAP)
    if Y_Arr(i) == 0 
        p1 = plot(cellAP(i).times,cellAP(i).V,'linewidth',0.3,'Color',c1);
        
    end
set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end
hold on
for i = 1:3:length(BaseCellsEAD)
        p2 = plot(BaseCellsEAD(i).times,BaseCellsEAD(i).V,'linewidth',0.3,'Color',c2);
        p12.Color(4) = 0.2;
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
settings.SubFolder = 'Dronedarone';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File,'CABaseCells');
load(File1,'Y_Arr')

arrhy = 1;
for i = 1:length(CABaseCells)
    ti = CABaseCells(i).times;
    Vi = CABaseCells(i).V;
    Caii = CABaseCells(i).Cai;
    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP(arrhy).times = t-1000;
    cellAP(arrhy).V = V;
    arrhy = arrhy + 1;
end
% 绘图
figure
for i = 1:length(cellAP)
    if Y_Arr(i) == 0 
        p1 = plot(cellAP(i).times,cellAP(i).V,'linewidth',0.3,'Color',c1);
        
    end
set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end
hold on
for i = 1:length(BaseCellsEAD)
        p2 = plot(BaseCellsEAD(i).times,BaseCellsEAD(i).V,'linewidth',0.3,'Color',c2);
        p12.Color(4) = 0.2;
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
settings.SubFolder = 'Ibutilide';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File,'CABaseCells');
load(File1,'Y_Arr')

arrhy = 1;
for i = 1:length(CABaseCells)
    ti = CABaseCells(i).times;
    Vi = CABaseCells(i).V;
    Caii = CABaseCells(i).Cai;
    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP(arrhy).times = t-1000;
    cellAP(arrhy).V = V;
    arrhy = arrhy + 1;
end
% 绘图
figure
for i = 1:length(cellAP)
    if Y_Arr(i) == 0 
        p1 = plot(cellAP(i).times,cellAP(i).V,'linewidth',0.3,'Color',c1);
        
    end
set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end
hold on
for i = 1:5:length(BaseCellsEAD)
        p2 = plot(BaseCellsEAD(i).times,BaseCellsEAD(i).V,'linewidth',0.3,'Color',c2);
        p12.Color(4) = 0.2;
hold on
end
hold on
for i = 1:7:length(BaseCellsEAD)
        p2 = plot(BaseCellsEAD(i).times,BaseCellsEAD(i).V,'linewidth',0.3,'Color',c2);
        p12.Color(4) = 0.2;
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
settings.SubFolder = 'Sotalol';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File,'CABaseCells');
load(File1,'Y_Arr')

arrhy = 1;
for i = 1:length(CABaseCells)
    ti = CABaseCells(i).times;
    Vi = CABaseCells(i).V;
    Caii = CABaseCells(i).Cai;
    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP(arrhy).times = t-1000;
    cellAP(arrhy).V = V;
    arrhy = arrhy + 1;
end
Y_Arr1 = Y_Arr;


load('TestPop/AFClassIII/male/Dronedarone/CABaseCells.mat','CABaseCells');
load('TestPop/AFClassIII/male/Dronedarone/Y_Arr.mat','Y_Arr')
arrhy = 1;
for i = 1:length(CABaseCells)
    ti = CABaseCells(i).times;
    Vi = CABaseCells(i).V;
    Caii = CABaseCells(i).Cai;
    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP2(arrhy).times = t-1000;
    cellAP2(arrhy).V = V;
    arrhy = arrhy + 1;
end
Y_Arr2 = Y_Arr;


% 绘图
figure
for i = 1:length(cellAP)
    if Y_Arr1(i) == 0 
        p1 = plot(cellAP(i).times,cellAP(i).V,'linewidth',0.3,'Color',c1);
        
    end
set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end

hold on
p2 = plot(BaseCellsEAD(7).times,BaseCellsEAD(7).V,'linewidth',0.3,'Color',c2);
p12.Color(4) = 0.2;
hold on
p2 = plot(BaseCellsEAD(19).times,BaseCellsEAD(19).V,'linewidth',0.3,'Color',c2);
p12.Color(4) = 0.2;
hold on
p2 = plot(BaseCellsEAD(31).times,BaseCellsEAD(31).V,'linewidth',0.3,'Color',c2);
p12.Color(4) = 0.2;
hold on
p2 = plot(BaseCellsEAD(29).times,BaseCellsEAD(29).V,'linewidth',0.3,'Color',c2);
p12.Color(4) = 0.2;
% 关闭右边和上边的坐标轴线
box off;
xlabel('Time (ms)','FontSize',6,'FontName','Calibri')
ylabel('Voltage (mV)','FontSize',6,'FontName','Calibri')
set(gcf,'Position',[219,161,448,392])
hold off
%}

%% Vernakalant
%
settings.SubFolder = 'Vernakalant';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File,'CABaseCells');
load(File1,'Y_Arr')

arrhy = 1;
for i = 1:length(CABaseCells)
    ti = CABaseCells(i).times;
    Vi = CABaseCells(i).V;
    Caii = CABaseCells(i).Cai;
    [t,V,~] = splitdata(Vi,Caii,ti);
    cellAP(arrhy).times = t-1000;
    cellAP(arrhy).V = V;
    arrhy = arrhy + 1;
end
% 绘图
figure
for i = 1:length(cellAP)
    if Y_Arr(i) == 0 
        p1 = plot(cellAP(i).times,cellAP(i).V,'linewidth',0.3,'Color',c1);
        
    end
set(gca,'FontSize',6,'FontName','Calibri','linewidth',0.5,...
        'YLim',[-100 50],'YTick',[-100 -50 0 50],'XLim',[0 1000],'XTick',[0 1000])
hold on
end
hold on
for i = 1:2:length(BaseCellsEAD)
        p2 = plot(BaseCellsEAD(i).times,BaseCellsEAD(i).V,'linewidth',0.3,'Color',c2);
        p12.Color(4) = 0.2;
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
