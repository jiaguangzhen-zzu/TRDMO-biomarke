function second_order_connectivity = calculate_second_order_connectivity(test_adjacency)
data=load('E:\数据\TCGA1\Data6_TCGA-CESC.mat');
for fnum=1
fieldName = ['N', num2str(fnum)];
        test_adjacency = data.(fieldName);
       
        Problem.D=size(test_adjacency,1); %维度
        test_net=zeros(Problem.D,Problem.D);
        test_net(test_adjacency~=0)=1;
end

%connectivity = sum(test_net, 2);
    num_genes = size(test_adjacency, 1);
    % 初始化一个数组，用于存储每个基因的二阶连通度
    second_order_connectivity = zeros(1, num_genes);
    
    % 计算每个基因的连接基因个数
    gene_connections = count_gene_connections(test_adjacency);
    sorted_onegene = sort(gene_connections, 'descend');
    % 遍历邻接矩阵的每一行或每一列

    for i = 1:num_genes
        % 找到第 i 行（或第 i 列）中非零元素的位置
        connected_indices = find(test_adjacency(:, i));
        % 计算与第 i 基因相连的基因的连接基因个数的总和
        total_connected_genes = gene_connections(i);
        % 统计这些相连基因相连的基因个数，并将它们加起来
        for j = 1:length(connected_indices)
            connected_gene = connected_indices(j);
            total_connected_genes = total_connected_genes + gene_connections(connected_gene);
        end
        % 将计算得到的二阶连通度赋值给相应的基因
        second_order_connectivity(i) = total_connected_genes;
    end
    n = size(test_adjacency, 1);  % 根据邻接矩阵的行数确定 n
    result_matrix = reshape(second_order_connectivity, n, 1);  % 将结果转换为 n*1 的矩阵

 

% 对 result_matrix 进行降序排序，并获取排序后元素的索引
[sorted_values, sorted_indices] = sort(result_matrix, 'descend');

% 输出排序后的结果和相应的索引
disp(sorted_values);  % 排序后的矩阵元素
disp(sorted_indices);  % 排序后的索引顺序



combined_matrix = horzcat(sorted_values, sorted_indices);



% 创建一个与您的矩阵相同大小的零矩阵
gene_matrix = zeros(5, num_genes);

% 获取前50%基因的索引
num_genes_to_set_one = ceil(0.5 * size(combined_matrix, 1)); % 前50%的行数
top_50_percent_indices = combined_matrix(1:num_genes_to_set_one, 2); % 获取第二列前50%行的索引





% 生成随机数矩阵，大小与 gene_matrix 相同
random_values(:, top_50_percent_indices) = rand(size(gene_matrix, 1), numel(top_50_percent_indices));

% 将随机数矩阵中大于0.2的位置设置为1
random_values(random_values > 0.3) = 1;
random_values(random_values <= 0.3) = 0;

gene_matrix = random_values;




end


function gene_connections = count_gene_connections(adjacency_matrix)
    % 获取邻接矩阵的阶数，即基因的数量
    num_genes = size(adjacency_matrix, 1);
    % 初始化一个数组，用于存储每个基因的连接基因个数
    gene_connections = zeros(1, num_genes);
    
    % 遍历邻接矩阵的每一行或每一列
    for i = 1:num_genes
        % 计算第 i 行（或第 i 列）中非零元素的个数
        gene_connections(i) = nnz(adjacency_matrix(:, i));
    end
end








