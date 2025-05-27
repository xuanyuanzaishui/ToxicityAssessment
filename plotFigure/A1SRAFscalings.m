clear

%% 输出有效药物的离子电导参数值，用于python绘制分布图
%同时输出均值和方差

%SR\male
load('TestPop\SR\male\SRpopICsm.mat') % population data
SRpopscalingsm(:,15) = 0; % Gender, male 0, female 1
SRpopscalingsm(:,16) = 0; % Class
popScalings = SRpopscalingsm;

mean1 = mean(SRpopscalingsm);
std1 = std(SRpopscalingsm);

%SR\female
load('TestPop\SR\female\SRpopICsf.mat') % population data
SRpopscalingsf(:,15) = 1; % Gender, male 0, female 1
SRpopscalingsf(:,16) = 0; % Class
popScalings = [popScalings;SRpopscalingsf];

mean2 = mean(SRpopscalingsf);
std2 = std(SRpopscalingsf);

%AF\male
load('TestPop\AF\male\AFpopICsm.mat') % population data
AFpopscalingsm(:,15) = 0; % Gender, male 0, female 1
AFpopscalingsm(:,16) = 1; % Class
popScalings = [popScalings;AFpopscalingsm];

mean3 = mean(AFpopscalingsm);
std3 = std(AFpopscalingsm);

%AF\female
load('TestPop\AF\female\AFpopICsf.mat') % population data
AFpopscalingsf(:,15) = 1; % Gender, male 0, female 1
AFpopscalingsf(:,16) = 1; % Class
popScalings = [popScalings;AFpopscalingsf];

mean4 = mean(AFpopscalingsf);
std4 = std(AFpopscalingsf);

meanStdMF = [mean1;std1;mean2;std2;mean3;std3;mean4;std4];

%% csv

outputlabels = {'GNa', 'GNaL', 'Gto', 'GKr', 'GKs', 'GKur',...
     'GK1','GCaL', 'Pnak', 'Gncx', 'GpCa', 'Rel', 'SERCA', 'Leak','Gender','class'};%labels

Output_Dir = 'E:\matlab\model\Table\'; % path

%SRAFscalingsMF
Output_File = 'SRAFscalingsMF.csv';
outputFile(Output_Dir,Output_File,popScalings,outputlabels)%特征

%SRAFscalingsMF
Output_File = 'SRAFscalingsMeanStdMF.csv';
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