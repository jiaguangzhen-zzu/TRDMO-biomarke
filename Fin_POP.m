path_in_RBM='E:\结果\Data23\FK\RBM';
path_in_RPRORBM='E:\结果\Data23\FK\RPRORBM';
path_in_MPMMEA='E:\对比试验结果\MPMMEA\Data23';
path_in_MOEA_PSL='E:\对比试验结果\MOEA_PSL\Data23';
path_in_DNSGA2='E:\对比试验结果\DNSGA2\Data23';

for fnum=1:4
        NDSol1=load([path_in_RBM,'\RBM' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']);
        NDSol2=load([path_in_RPRORBM,'\RPRORBM0.2四' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']);
        NDSol3=load([path_in_MPMMEA,'\MPMMEA' num2str(fnum) '_MPMMEA-RBM_NDSol.mat']);
        NDSol4=load([path_in_MOEA_PSL,'\MOEA_PSL' num2str(fnum) '_MPMMEA-RBM_NDSol.mat']);
        NDSol5=load([path_in_DNSGA2,'\DNSGA2一' num2str(fnum) '_MMPDENB-RBM_NDSol.mat']);

% 假设有一个 cell 数组存储多个结构体
NDsol = {NDSol1, NDSol2, NDSol3, NDSol4, NDSol5};  % 将你的结构体放在一个 cell 数组中

% 设置颜色
colors = {'r', 'g', 'b', 'o', 'm'}; % 红色, 绿色, 蓝色, 橙色, 紫色

% 创建新的图形
figure;
hold on;

for i = 1:length(NDsol)
    matrix = NDsol{i}.FVV;  % 获取 FVV 矩阵
    
    % 处理 y 轴坐标
    y_processed = log2(matrix(:, 2) + 2);  % 对 y 加 2 后取 log2
    y_processed = 1 ./ y_processed;          % 对 y 取倒数
    
    % 处理 x 轴坐标
    x_processed = matrix(:, 1);               % 直接获取 x 坐标
    valid_indices = x_processed <= 100;       % 找到有效的 x 坐标
    
    % 过滤 x 和 y
    x_filtered = x_processed(valid_indices);
    y_filtered = y_processed(valid_indices);
    
    % 对 x 和 y 进行归一化
    x_normalized = (x_filtered - min(x_filtered)) / (max(x_filtered) - min(x_filtered));
    y_normalized = (y_filtered - min(y_filtered)) / (max(y_filtered) - min(y_filtered));
    
    % 绘制处理后的点
    scatter(x_normalized, y_normalized, 100, colors{i}, 'filled', 'DisplayName', ['Cell ' num2str(i)]);
end

% 添加标签和图例
xlabel('Normalized X-axis');
ylabel('Normalized Y-axis');
title('Processed Points from NDsol Cells');
legend;
grid on;
hold off;
        
end