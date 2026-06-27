function score=Cal_HV_my()
%data=load('E:\数据\TCGA1\Data10_TCGA-ESCA.mat');
base_path = 'E:\结果\Data11\FK\比较单次';
folder_template = 'Folder%d';
results_cell = cell(30, 4);
for folder_index = 1:30
    % 构建当前文件夹的完整路径
    path_in_TRPRORBM = fullfile(base_path, 'TRPRORBM', '30000', sprintf(folder_template, folder_index));
    path_in_RPRORBM = fullfile(base_path, 'RPRORBM', '30000', sprintf(folder_template, folder_index));
    path_in_TPRORBM = fullfile(base_path, 'TPRORBM', '30000', sprintf(folder_template, folder_index));
    path_in_MMPDENBRBM = fullfile(base_path, 'RBM', '30000', sprintf(folder_template, folder_index));
    path_in_PRORBM = fullfile(base_path, 'PRORBM', '30000', sprintf(folder_template, folder_index));
     for fnum = 1 : 4
          fieldName = ['N', num2str(fnum)];
        pop1_data = load(fullfile(path_in_TRPRORBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_NDSol.mat']));
        pop2_data = load(fullfile(path_in_RPRORBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_NDSol.mat']));
        pop3_data = load(fullfile(path_in_TPRORBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_NDSol.mat']));
        pop4_data = load(fullfile(path_in_MMPDENBRBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_NDSol.mat']));
        pop5_data = load(fullfile(path_in_PRORBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_NDSol.mat']));    
        results_cell{folder_index, fnum} = {pop1_data.FVV, pop2_data.FVV, pop3_data.FVV, pop4_data.FVV, pop5_data.FVV};
     end
end
for col = 1:4
    column_data = results_cell(:, col);
    result = cell(1, 5);

    for col = 1:5
        current_column_data = cell(30, 1);
        for row = 1:30
            current_column_data{row} = column_data{row}(col);
        end
        result{col} = current_column_data;
    end

    colors = ['r', 'g', 'b', 'm', 'c']; % 红、绿、蓝、品红、青
    figure;

    allX = [];
    allY = [];
    
    % 收集所有数据
    for col = 1:5
        matrices = result{col};
        for i = 1:length(matrices)
            data = matrices{i};
            data = data{1};
            allX = [allX; data(:, 1)];
            allY = [allY; data(:, 2)];
        end
    end

    % 计算异常值阈值（前5%）
    thresholdIndexX = ceil(0.05 * length(allX));
    thresholdIndexY = ceil(0.05 * length(allY));
    
    sortedX = sort(allX);
    sortedY = sort(allY);

    cutoffX = sortedX(thresholdIndexX); % X轴异常值阈值
    cutoffY = sortedY(thresholdIndexY); % Y轴异常值阈值

    % 筛选有效数据（去除异常值）
    validIndices = (allX > cutoffX) & (allY > cutoffY); 
    filteredX = allX(validIndices);
    filteredY = allY(validIndices);

    % 对 Y 轴数据进行 log2 变换，只保留大于 0 的值
    logY = log2(filteredY);
    validLogY = logY(logY > 0); % 只保留大于 0 的值
    
    % 对应的 X 轴数据
    validFilteredX = filteredX(logY > 0); % 只保留对应的 X 数据

    % 取倒数
    invLogY = 1 ./ validLogY;

    % 归一化 X 和 Y 数据
    maxX = max(validFilteredX);
    minX = min(validFilteredX);
    maxY = max(invLogY);
    minY = min(invLogY);

    % 绘制归一化图形
    for col = 1:5
        matrices = result{col};
        for i = 1:length(matrices)
            data = matrices{i};
            data = data{1};
            currentX = data(:, 1);
            currentY = data(:, 2);
            
            % 仅保留有效数据
            validIndices = (currentX > cutoffX) & (currentY > cutoffY);
            validX = currentX(validIndices);
            validY = currentY(validIndices);

            % 对 Y 进行 log2 变换和倒数处理
            logY = log2(validY);
            validLogY = logY(logY > 0);
            validFilteredX = validX(logY > 0); % 只保留对应的 X 数据
            invLogY = 1 ./ validLogY;

            % 对 X 和 Y 进行归一化
            normX = (validFilteredX - minX) / (maxX - minX);
            normY = (invLogY - min(min(invLogY))) / (maxY - min(min(invLogY)));

            plot(normX, normY, 'o', 'Color', colors(col), ...
                'DisplayName', ['Group ', num2str(col)]);
            hold on;
        end
    end

    legend show;
    xlabel('归一化节点');
    ylabel('归一化分数（经过log2和取倒数处理）');
    title('PF');
    hold off; 
end
% for col = 1:4
%     column_data = results_cell(:, col);
%     result = cell(1, 5);
% 
%     for col = 1:5
%         current_column_data = cell(30, 1);
%         for row = 1:30
%             current_column_data{row} = column_data{row}(col);
%         end
%         result{col} = current_column_data;
%     end
%    
% 
%     
%     colors = ['r', 'g', 'b', 'm', 'c']; % 红、绿、蓝、品红、青
%     figure;
%     for col = 1:5
%         matrices = result{col};
%         for i = 1:length(matrices)
%             data = matrices{i};
%             data = data{1};
%             plot(data(:, 1), data(:, 2), 'o', 'Color', colors(col), ...
%                 'DisplayName', ['Group ', num2str(col)]);
%             hold on;
%         end
%     end
%    legend show;
%     xlabel('节点');
%     ylabel('分数');
%     title('PF');
%     hold off; 
% 
% end