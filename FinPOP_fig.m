% data=load('E:\数据\TCGA降维矩阵\Data27_TCGA-THCA.mat');
% path_in_TRDMO='E:\结果\Data27\FK\RPRORBM降维';
% path_in_RBM='E:\结果\Data27\FK\RBM降维';
% path_in_MOEA_PSL='E:\对比试验结果\MOEA_PSL\Data27';
% path_in_MPMMEA='E:\对比试验结果\MPMMEA\Data27';
% path_in_DNSGA2='E:\对比试验结果\DNSGA2\Data27';
% path_in_MP_DFR='E:\对比试验结果\MP_DFR\Data27';
% path_in_HATC='E:\对比试验结果\HATC\Data27';
%    for fnum =  3:4 %length(numf)
% fieldName = ['N', num2str(fnum)];
%         test_adjacency = data.(fieldName);
%         NDSol1=load([path_in_TRDMO,'\RPRORBM' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']);
%         NDSol1=NDSol1.FVV;
%         NDSol2=load([path_in_RBM,'\RBM' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']);
%         NDSol2=NDSol2.FVV;
%         NDSol3=load([path_in_MOEA_PSL ,'\MOEA_PSL降维' num2str(fnum) '_MPMMEA-RBM_NDSol.mat']);
%         NDSol3=NDSol3.FVV;
%         NDSol4=load([path_in_MPMMEA,'\MPMMEA降维' num2str(fnum) '_MPMMEA-RBM_NDSol.mat']);
%         NDSol4=NDSol4.FVV;
%         NDSol5=load([path_in_DNSGA2,'\DNSGA2降维' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']);
%         NDSol5=NDSol5.FVV;
%         NDSol6=load([path_in_MP_DFR,'\降维' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']);
%         NDSol6=NDSol6.FVV;
%         NDSol7=load([path_in_HATC,'\降维' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']); 
%         NDSol7=NDSol7.FVV;
%    end


% 基础配置（关闭TeX，_正常显示）
set(0,'DefaultAxesFontName','SimHei');
set(0,'DefaultTextFontName','SimHei');
set(0,'DefaultLegendFontName','SimHei');
set(0,'DefaultTextInterpreter','none');

% 1. 指定fnum
fnum = 3;

% 2. 加载数据
NDSol1 = load(['E:\结果\Data7\FK\RPRORBM降维\RPRORBM3',num2str(fnum),'_MMPDENB-RBM_NDSol.mat']).FVV;
NDSol2 = load(['E:\结果\Data7\FK\RBM降维\RBM',num2str(fnum),'_MMPDENB-RBM_NDSol.mat']).FVV;
NDSol3 = load(['E:\对比试验结果\MOEA_PSL\Data7\MOEA_PSL降维',num2str(fnum),'_MPMMEA-RBM_NDSol.mat']).FVV;
NDSol4 = load(['E:\对比试验结果\MPMMEA\Data7\MPMMEA降维',num2str(fnum),'_MPMMEA-RBM_NDSol.mat']).FVV;
NDSol5 = load(['E:\对比试验结果\DNSGA2\Data7\DNSGA2降维',num2str(fnum),'_MMPDENB-RBM_NDSol.mat']).FVV;
NDSol6 = load(['E:\对比试验结果\MP_DFR\Data7\降维',num2str(fnum),'_MMPDENB-RBM_NDSol.mat']).FVV;
NDSol7 = load(['E:\对比试验结果\HATC\Data7\降维',num2str(fnum),'_MMPDENB-RBM_NDSol.mat']).FVV;

% 3. 定义颜色和名称
colors = [1 0 0; 1 0.6 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 0.5 0 0.5];
names = {'TRDMO','MPDENB-RBM','MOEA\_PSL','MPMMEA','DNSGA-II','MP\_DFR','HATC'};
allData = {NDSol1,NDSol2,NDSol3,NDSol4,NDSol5,NDSol6,NDSol7};

% 4. 绘图（点大小设为50）
figure('Position',[100,100,900,600]); hold on;
for i=1:7
    if ~isempty(allData{i})

 tempData = allData{i};
tempData = tempData(tempData(:, 1) <= 100, :); % 核心：删除第一列>100的整行

% 第二步：从筛选后的矩阵取列（此时col1已无>100的元素）
col1 = tempData(:,1);
col2_raw = tempData(:,2);
col2_processed = zeros(size(col2_raw));
    
        for j=1:length(col2_raw)
            if col2_raw(j) > 2
                col2_processed(j) = 1 / (log2(col2_raw(j)) / 10);
            else
                col2_processed(j) = 1 / (col2_raw(j) / 10);
            end
        end
        
        scatter(col1,col2_processed,80,colors(i,:),'filled','MarkerFaceAlpha',0.6);
    end
end

set(gca,'YScale','log');
xlabel('第一目标值（已最小化）');
ylabel('处理后第二目标值（越小越好）');
title(['fnum=',num2str(fnum),' 7方法目标值散点图']);

% 核心：调整图例大小（文字12号，标记10号，整体缩放1.1）
lgd = legend(names,'Location','best','FontSize',15);
lgd.EntryContainer.Children.MarkerSize = 80; % 图例内点标记大小
lgd.Scale = 1.1; % 图例整体放大10%

grid on; grid alpha 0.3;
hold off;