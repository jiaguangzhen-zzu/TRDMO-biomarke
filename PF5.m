function score=Cal_HV_my2()
%% 701 4种算法对同一个病人使用相同的参考点，使用秩和检验来判断好坏 length(numf)
%% 参数敏感性分析，对随机抽取的病人计算HV值
data=load('E:\数据\TCGA1\Data10_TCGA-ESCA.mat');
base_path = 'E:\结果\Data10\FK\比较单次';
folder_template = 'Folder%d';
results_cell = cell(30, 4);
% score1=zeros(4,1);
% score2=zeros(4,1);
% 
% std1=zeros(4,1);
% std2=zeros(4,1);

for folder_index = 1:30
    path_in_My = fullfile(base_path, 'TRPRORBM', '30000', sprintf(folder_template, folder_index));
    path_in_RPRORBM = fullfile(base_path, 'RPRORBM', '30000', sprintf(folder_template, folder_index));
    path_in_TPRORBM = fullfile(base_path, 'TPRORBM', '30000', sprintf(folder_template, folder_index));
    path_in_MMPDENBRBM = fullfile(base_path, 'RBM', '30000', sprintf(folder_template, folder_index));
    path_in_PRORBM = fullfile(base_path, 'PRORBM', '30000', sprintf(folder_template, folder_index));
    for fnum = 1 : 4%length(numf)
        fieldName = ['N', num2str(fnum)];
        test_adjacency = data.(fieldName);
        %% 获取参考点
        pop1_data = load(fullfile(path_in_MMPDENBRBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_POP.mat']));    
        pop2_data = load(fullfile(path_in_PRORBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_POP.mat']));
        pop3_data = load(fullfile(path_in_RPRORBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_POP.mat']));
        pop4_data = load(fullfile(path_in_TPRORBM, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_POP.mat']));
        pop5_data = load(fullfile(path_in_My, ['NEGRBM', num2str(fnum), '_MMPDENBFK-RBM_POP.mat']));
        results_cell{folder_index, fnum} = {pop1_data.Final_Pop, pop2_data.Final_Pop, pop3_data.Final_Pop, pop4_data.Final_Pop, pop5_data.Final_Pop};
       
    end

end
for l = 1:4
    fieldName = ['N', num2str(l)];
    test_adjacency1 = data.(fieldName);
    column_cells = results_cell(:, l);
    merged_cell{l} = vertcat(column_cells{:});
    merged_cell1 =merged_cell(:, l);
    colors = ['r', 'g', 'b', 'm', 'c']; % 红、绿、蓝、品红、青
    for col = 1:5
        POP = merged_cell1{1, 1} (:, col);
        total_rows = sum(cellfun(@(x) size(x, 1), POP));
        num_cols = size(POP{1}, 2);
        pop = zeros(total_rows, num_cols);
        current_row = 1;
        for i = 1:30
            matrix = POP{i};
            num_rows = size(matrix, 1);
            pop(current_row:(current_row + num_rows - 1), :) = matrix;
            current_row = current_row + num_rows;
        end
       [pop,~,~]=unique(pop,'rows');
       functionvalue = Calfunctionvalue(pop,test_adjacency1);
       [FrontNo,~] = NDSort(functionvalue,size(functionvalue,1));      
       FV=functionvalue(FrontNo==1,:);
       FV(:, 2) = -FV(:, 2);
       FV(:, 2) = 1 ./ FV(:, 2);
       plot(FV(:,1), FV(:,2), 'o', 'Color', colors(col), 'DisplayName', ['Column ' num2str(col)]);
      hold on; % 保持图像，继续在同一图上绘制
    end
    legend('show');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Visualization of FV for Each Column');
    hold off; % 释放图像


end
end

function f=calhv(population,Ref_point,test_adjacency)
population = population.Final_Pop;


population11 = population;
f1 = Calfunctionvalue(population11,test_adjacency);
% f1(:,2)=f1(:,2)+add;
f1(:,2)=f1(:,2);
f=HV(f1,Ref_point);

if f>1
    f=0;
end
%end
end

%% 
function Functionvalue = Calfunctionvalue(pop,test_adjacency)
%pop = pop.Final_Pop;
Functionvalue=zeros(size(pop,1),2);
functionvalue=zeros(size(pop,1),4);
for i=1:size(pop,1)
    Functionvalue(i,1)=sum(pop(i,:));
    %计算内部矩阵
    a=find(pop(i,:)==1);
    matrix=test_adjacency(a,:);
    inmatrix=matrix(:,a);
    genin=nonzeros(inmatrix);

    functionvalue(i,1)=abs(mean(genin));
    functionvalue(i,2)=std(genin);

    %计算外部矩阵
    gen=nonzeros(matrix);
    genout=setdiff(gen,genin);

    functionvalue(i,3)=abs(mean(genout));
    functionvalue(i,4)=functionvalue(i,1)*functionvalue(i,2)/functionvalue(i,3);

end
Functionvalue(:,2)= -round(functionvalue(:,4)*100)/100;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,:)=[];
ttt= Functionvalue(:,2)==0;
Functionvalue(ttt,:)=[];
end
%% 




