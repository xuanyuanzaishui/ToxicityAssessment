clear
%% 该函数用来输出SR AF的特征csv，用于python绘制分布图
%同时输出均值和方差
load('TestPop\SR\male\SRfeaturesm.mat') % SR male population feture
load('TestPop\SR\female\SRfeaturesf.mat') % SR female population feture

%load('TestPop\AF\male\AFfeaturesm.mat') % AF male population feture
%load('TestPop\AF\female\AFfeaturesf.mat') % AF female population feture

% SR male
mean1 = mean(SRfeaturesm);
std1 = std(SRfeaturesm);
SRfeaturesm(:,15) = 0; % gender
SRfeaturesm(:,16) = 0; % class

features = SRfeaturesm;

% SR female
mean2 = mean(SRfeaturesf);
std2 = std(SRfeaturesf);
SRfeaturesf(:,15) = 1; % gender
SRfeaturesf(:,16) = 0; % class

features = [features;SRfeaturesf];

meanStdMFsr = [mean1;std1;mean2;std2;];

% AF male
mean3 = mean(AFfeaturesm);
std3 = std(AFfeaturesm);
AFfeaturesm(:,15) = 0; % gender
AFfeaturesm(:,16) = 1; % class

features = [features;AFfeaturesm];

% AF female
mean4 = mean(AFfeaturesf);
std4 = std(AFfeaturesf);
AFfeaturesf(:,15) = 1; % gender
AFfeaturesf(:,16) = 1; % class

features = [features;AFfeaturesf];
meanStdMF = [mean1;std1;mean2;std2;mean3;std3;mean4;std4];

%% csv

outputlabels = {'RMP', 'dV/dtmax', 'APA', 'APD20', 'APD40', 'APD50',...
     'APD90','APDtri', 'CTA', 'CTmax', 'CTD50', 'CTD90', 'CTDtri', 'CTD','Gender','class'};%labels

Output_Dir = 'E:\matlab\model\Table\'; % path

%SRAFfeatureMF
Output_File = 'SRAFfeatureMF.csv';
outputFile(Output_Dir,Output_File,features,outputlabels)%特征

%SRAFmeanstdMF
Output_File = 'SRAFfeantureMeanStdMF.csv';
outputFile(Output_Dir,Output_File,meanStdMF,outputlabels)%特征



%% hanshu

%输出文件函数
function outputFile(Output_Dir,Output_File,outputFeature,outputlabels)
        path = fullfile(Output_Dir,Output_File);
        if isfile(path)
            disp('FILE ALREADY EXISTS.')
        else
            fid = fopen(Output_File, 'w') ;
            fprintf(fid, '%s,', outputlabels{1,1:end-1}) ;
            fprintf(fid, '%s\n', outputlabels{1,end}) ;
            fclose(fid) ;
            dlmwrite(Output_File, outputFeature, '-append') ;
        end
end
%}