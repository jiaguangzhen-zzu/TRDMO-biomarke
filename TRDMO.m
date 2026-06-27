function [outputArg1,outputArg2] = MMPDENB_PSL(inputArg1,inputArg2)
%% MOEA-PSL identify edge-network biomarkers
%   此处显示详细说明
data=load('E:\数据\TCGA1\Data24_TCGA-SKCM.mat');

    experiment_num=1;
    Problem.maxFE=30;
    Problem.N=300;
    Problem.lower=0;
    Problem.upper=1;

 population_cell = cell(4, 30);
 for EXP_NUM=1:experiment_num
     current_populations = cell(1, 4);
    for fnum=1:4 %length(Sample_net)%邻接矩阵
        fieldName = ['N', num2str(fnum)];
        test_adjacency = data.(fieldName);       
        Problem.D=size(test_adjacency,1); %维度
        
        test_net=zeros(Problem.D,Problem.D);
        test_net(test_adjacency~=0)=1;
        %% paired edge :g
        Cons=CONS(test_net);
        Cons(all(Cons==0,2),:)=[];      
        Cnum1=size(Cons,1);
        g=cell(Cnum1,1);     
        for i=1:Cnum1
            g{i}=find(Cons(i,:));          
        end     
        Cnum=size(g,1);
        %clear Cons
        %% 计算每个维度的目标函数值
        Dimension=eye(Problem.D,Problem.D);
        CV_num=Calcons(Problem.D,Cnum,g,Dimension);
        D_score=abs(CV_num-Cnum);
        %Non_dominated_sol=cell(30,1);
        R=cell(30,1);
        %for EXP_NUM=1:experiment_num
        CalNum=0;
        max_value = max(D_score);
        scaled_value = max_value * 0.2;
        rounded_value = round(scaled_value);
       top_indices = find( D_score > rounded_value);    
            %% Population initialization
            REAL = false;

            if  fnum == 1
                Dec = ones(Problem.N,Problem.D);
                Mask=creatpop(Problem.N,Problem.D,D_score,g);
                [Mask,~] = unique(Mask,'rows');
                if size(Mask,1)<Problem.N
                    newsoultion=Problem.N-size(Mask,1);
                    ns=creatpop(newsoultion,Problem.D,D_score,g);
                    Mask=[Mask;ns];
                end

                Population = Dec.*Mask;
                [Population,calnum]  = SOLUTION(Population,test_adjacency,g,D_score);
                CalNum=CalNum+calnum;

                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);%%里面的非支配排序要换成自己的
            else
                Dec = ones(Problem.N,Problem.D);


                if size( SelectedPop, 1) < 100

                    neededRows = 100 - size(SelectedPop, 1);
                    part1 = SelectedPop;
                    part2 = creatpop(neededRows,Problem.D,D_score,g); % 生成补足的随机数据
                else

                    part1 =SelectedPop(1:100, :);
                    part2 = [];  % 不需要额外的补足行
                end
                part3 = creatpop(200,Problem.D,D_score,g);

                part = [part1; part2];
                [Part3,~]  = SOLUTION(part3,test_adjacency,g,D_score);
                %CalNum=CalNum+calnum;
                x_values = zeros(1, 200);
                x_values_jd = zeros(1, 200);
                for i = 1:200
                    obj_values = Part3{1, i}.obj(:, 2);
                    x_values(i) = obj_values;
                    x_jd= Part3{1, i}.obj(:, 1);
                    x_values_jd(i) = x_jd;
                end
                [Part,~]  = SOLUTION(part,test_adjacency,g,D_score);
                %CalNum=CalNum+calnum;
                y_values = zeros(1, 100);
                y_values_jd = zeros(1, 100);
                for i = 1:100
                    obj_values = Part{1, i}.obj(:, 2);
                    y_values(i) = obj_values;
                    y_jd= Part{1, i}.obj(:, 1);
                    y_values_jd(i) = y_jd;
                end
                max1 = max(abs(x_values));
                max2 = max(abs(y_values));
                max_all = max(max1, max2);
                fvv1 = abs(x_values) / max_all;
                fvv2 = abs(y_values )/ max_all;
                max3 = max(1 ./x_values_jd);
                max4 = max(1 ./y_values_jd);
                max_all1 = max(max3, max4);
                jd1 = abs(1 ./x_values_jd) / max_all1;
                jd2 = abs(1 ./y_values_jd )/ max_all1;
                lastpart = cell(1, 9);
                for g = 1:9
                    weight1 = 0.1 * g;
                    weight2 = 1 - weight1;
                    %sums1(i, :) = weight1 * fvv1 + weight2 * jd1;
                    %sums2(i, :) = weight1 * fvv2 + weight2 * jd2;
                    sums1 = weight1 * fvv1 + weight2 * jd1;
                    sums2 = weight1 * fvv2 + weight2 * jd2;
                    averagesum = mean(sums1);
                    repair_indices = sums2 < averagesum;
                    stay_indices = sums2 >= averagesum;
                    stay_populations = part(stay_indices, :);
                    repair_populations = part(repair_indices, :);
                    %random_indicate =randperm(100,50);
                    %random_population = part(random_indicate,:);
                 random_population = part(repair_indices,:);
                    for p = 1:size(random_population)
                        %p = 1:length( repair_indices)
                        population_data = random_population(p, :);
                        original_ones = population_data == 1;
                        new_population = zeros(size(population_data));
                        new_population(original_ones) = 1;
                        for i = 1:length(top_indices)
                            index = top_indices(i);
                            if population_data(index) == 0
                                if rand <= 0.5
                                    new_population(index) = 1;
                                end
                            else
                                new_population(index) = 1;
                            end
                            repair_populations(p, :) = new_population;
                        end

                    end
                    lastpart{g} = [ stay_populations;repair_populations];
                end
                combined_matrix = vertcat(lastpart{:});
                [combined_matrix,~,~]=unique(combined_matrix,'rows');
                functionvalue = Calfunctionvalue(combined_matrix,test_adjacency);
                [FrontNo,~] = NDSort(functionvalue,size(functionvalue,1));
                [~, ParetoFrontRank] = sort(FrontNo);
                ParetoFrontIndices = ParetoFrontRank(1:100);
                repairedpop = combined_matrix(ParetoFrontIndices, :);
                Mask = [repairedpop; part3];  % 总共300行
                if size(Mask,1)<Problem.N
                    newsoultion=Problem.N-size(Mask,1);
                    ns=creatpop(newsoultion,Problem.D,D_score,g);
                    Mask=[Mask;ns];
                end
                Population = Dec.*Mask;
                [Population,calnum]  = SOLUTION(Population,test_adjacency,g,D_score);
                CalNum=CalNum+calnum;
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);%%里面的非支配排序要换成自己的
            end
            %% Optimization
            rho = 0.5;
            RHO=0.5;
            while(CalNum<Problem.maxFE)
                if CalNum<ceil(Problem.maxFE/2) %%前半段用原始空间去尽可能的收敛
                    Site=false(1,ceil(Problem.N));
                    if any(Site)
                        [rbm,~,allZero,allOne] = ModelTraining(Mask(index,:),Dec(index,:),REAL);%% 这才是用NP解去训练RBM
                    else
                        [rbm,~,allZero,allOne] = deal([]);
                    end
                    
                    Popmarks   = extract_d(Population);
                    MatingPool = TournamentSelection_hamming(length(Population),length(Population),Popmarks,FrontNo,CrowdDis);
                    [OffMask,poss_s_num] = Operator1(MatingPool,rbm,Site,allZero,allOne,D_score,test_adjacency);
                    OffMask =logical(OffMask);
                    OffDec  =ones(size(OffMask,1),size(OffMask,2));
                    Offspring = OffDec.*OffMask;
                    [Offspring,calnum]  = SOLUTION(Offspring,test_adjacency,g,D_score);
                    CalNum=CalNum+calnum;
                    
                    [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),poss_s_num);
                    rho = (rho+sRatio)/2;
                    RHO=[RHO;rho];
                else %% 后半段在充分了解pf的情况下，用RBM去学习pareto最优子空间
                    if rho<0.5
%                         rho=0.5;
                        index=FrontNo<ceil(max(FrontNo)/2);
                    else
                        index=FrontNo==1;
                    end
                    Site = rho > rand(1,ceil(Problem.N));  %% 只是表明多少个后代要用RBM产生
                    
                    if any(Site)
                        [rbm,~,allZero,allOne] = ModelTraining(Mask(index,:),Dec(index,:),REAL);%% 这才是用NP解去训练RBM

                    else
                        [rbm,~,allZero,allOne] = deal([]);
                    end
                    
                    Popmarks   = extract_d(Population);
                    MatingPool = TournamentSelection_hamming(length(Population),length(Population),Popmarks,FrontNo,CrowdDis);
                    [OffMask,poss_s_num] = Operator1(MatingPool,rbm,Site,allZero,allOne,D_score,test_adjacency);
                    OffMask =logical(OffMask);
                    OffDec  =ones(size(OffMask,1),size(OffMask,2));
                    Offspring = OffDec.*OffMask;
                    [Offspring,calnum]  = SOLUTION(Offspring,test_adjacency,g,D_score);
                    CalNum=CalNum+calnum;
                    
                    [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),poss_s_num);
                    rho = (rho+sRatio)/2;
                    RHO=[RHO;rho];
                end
            end
            PopObj      = extract(Population);
            Pop         = extract_d(Population);
            [FrontNo,~] = NDSort(PopObj,size(PopObj,1));
            outputpop = Pop(FrontNo==1,:);   

               %ParetoFront = 1;  
                  %ParetoFrontIndices = find(FrontNo == ParetoFront); 
                 %if length(ParetoFrontIndices) > 100
                    %ParetoFrontIndices = ParetoFrontIndices(1:100);  
                 %end
                    %SelectedPop = Pop(ParetoFrontIndices, :);
           % 计算每个个体的 Pareto 前沿等级

           [~, ParetoFrontRank] = sort(FrontNo);
           ParetoFrontIndices = ParetoFrontRank(1:100);
           SelectedPop = Pop(ParetoFrontIndices, :);

           current_populations{fnum} = outputpop();
       
     end
      for fnum = 1:4
        population_cell{fnum,EXP_NUM} = current_populations{fnum};
      end

  end
 %end
  Population_cell = cell(1, 4);
       for row = 1:4
   
    matrices_in_row = population_cell(row, :);
    concatenated_matrix = vertcat(matrices_in_row{:});
   Population_cell{row} = concatenated_matrix;
       end


     transposed_cell = population_cell' ;  
     for fnum = 1:4
        fieldName = ['N', num2str(fnum)];
        test_adjacency1 = data.(fieldName);   
        Non_dominated_sol1=Population_cell{fnum};
        Non_dominated_sol = transposed_cell(:, fnum);

      
        [pop,~,~]=unique(Non_dominated_sol1,'rows');%将得到的非支配解进行去重
        
        [functionvalue,~] = Calfunctionvalue(pop,test_adjacency1);%计算去重后的非支配解的目标函数
              
         
        
        [FrontNo,~] = NDSort(functionvalue,size(functionvalue,1));%对这些解进行pareto等级排序
        
        POP=pop(FrontNo==1,:);%输出最终的非支配解对应个体
        
        FV=functionvalue(FrontNo==1,:);
        
        [FVV,mod_position]=sortrows(FV);
        
        Final_Pop=POP(mod_position,:);% 输出最终排序过的个体
     

      
        %% 统计多模态出现的次数
        FVV(:,2)=-FVV(:,2);
        A=tabulate(FVV(:,1));
        B= A(:,2)~=0;
        c=A(B,:);
        %% 储存数据
            
        
        filename=strcat('E:\结果\Data24\FK\RPRORBM\RPRORBM0.2随机修补一半迁移种群',num2str(fnum),'_MMPDENB-RBM_NDSol.mat');
        save(filename,'FVV');
        filename1=strcat('E:\结果\Data24\FK\RPRORBM\RPRORBM0.2随机修补一半迁移种群',num2str(fnum),'_MMPDENB-RBM_Num.mat');
        save(filename1,'c');
        filename2=strcat('E:\结果\Data24\FK\RPRORBM\RPRORBM0.2随机修补一半迁移种群',num2str(fnum),'_MMPDENB-RBM_POP.mat');
        save(filename2,'Final_Pop');
        filename3=strcat('E:\结果\Data24\FK\RPRORBM\RPRORBM0.2随机修补一半迁移种群',num2str(fnum),'_MMPDENB-RBM_boxchart.mat');
        save(filename3,'Non_dominated_sol');
        filename4=strcat('E:\结果\Data24\FK\RPRORBM\RPRORBM0.2随机修补一半迁移种群',num2str(fnum),'_MMPDENB-RBM_pho.mat');
        save(filename4,'R');
    
    end
  %end
end


function pop=creatpop(popnum,D,D_score,g)
%% 选择两个一定相连的基因
pop=zeros(popnum,D);
gg=cell2mat(g);
%% 剩下的随机选择
for i=1:popnum
    mustselectnum_position=randperm(size(gg,1),1);  %% 选择两个一定相连的基因置为1
    mustselectnum=gg(mustselectnum_position,:);
    pop(i,mustselectnum)=1;
    D_score(mustselectnum)=0;
    canditateD=find(D_score~=0);
    
    min_D_score=D_score;
    gennum=round(0.9*rand*size(canditateD,1));
    for j=1:gennum
        canditateD=find(min_D_score~=0);   %% 使之前选过的基因不在参与选择
        variables=randperm(size(canditateD,1),2);
        if D_score(canditateD(variables(1)))>D_score(canditateD(variables(2)))
            pop(i,canditateD(variables(1)))=1;
            
            min_D_score(canditateD(variables(1)))=0;
        else
            pop(i,canditateD(variables(2)))=1;
            
            min_D_score(canditateD(variables(2)))=0;
        end
    end
end
end
function A_adjacent=CONS(test_Net)
[z1,z2]=find(triu(test_Net)~=0);
z=[z1,z2];
NNN=length(test_Net);

N1=NNN;
[N2,~]=size(z);
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N2,N1);
for i=1:N2
    
    A_adjacent(i,z(i,1))=1;
    A_adjacent(i,z(i,2))=1;
    
end
end
function CVvalue = Calcons(popnum,Cnum,g,pop)
cv=zeros(Cnum,1);
CVvalue=zeros(popnum,1);
for i=1:popnum
    ind=pop(i,:);
    for j=1:Cnum
        cv(j)=max(1-sum(ind(g{j})),0);
    end
    CVvalue(i)=sum(cv);
end
end
function [Population,calnum]  = SOLUTION(Mask,test_adjacency,g,D_score)
[functionvalue,calnum1] = Calfunctionvalue(Mask,test_adjacency);
[Mask,tt]=FIX2(Mask,functionvalue,g,size(test_adjacency,2),D_score,test_adjacency);
[functionvalue,calnum2]=Calfunctionvalue_afterfix(Mask,test_adjacency,tt,functionvalue);
calnum=calnum1+calnum2;
Population=cell(1,size(Mask,1));
for i=1:size(Population,2)
    Population{1,i}.dec=Mask(i,:);
    Population{1,i}.obj=functionvalue(i,:);
    Population{1,i}.con=0;
    Population{1,i}.add=[];
end
end
function [Functionvalue, calnum]= Calfunctionvalue(pop,test_adjacency)
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
calnum=i;
Functionvalue(:,2)= -round(functionvalue(:,4)*100)/100;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,1)=size(test_adjacency,2);
Functionvalue(tt,2)=0;
end
function [Functionvalue,calnum] = Calfunctionvalue_afterfix(pop,test_adjacency,tt,functionvalue)

for i=1:size(tt,1)
    functionvalue(tt(i),:)=Calfunctionvalue(pop(tt(i),:),test_adjacency);
end
calnum=i;
Functionvalue=functionvalue;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,2)=0;
ttt= Functionvalue(:,2)==0;
Functionvalue(ttt,1)=size(test_adjacency,2);
end
function Parents = TournamentSelection_hamming(K,N,Population,FrontNo,SpCrowdDis)
index=zeros(K,1);
% Parents=zeros(K,size(Population,2));
hamming_dist=pdist2(Population,Population,'hamming');
[~,site2]=sort(hamming_dist,2);
K1=randperm(N,K);  %第一个父母编号
K2=site2(K1,2);    %第二个父母编号

for i=1:K
    if FrontNo(K1(i))<FrontNo(K2(i))
        index(i)=K1(i);
    elseif FrontNo(K1(i))==FrontNo(K2(i))
        if SpCrowdDis(K1(i))>=SpCrowdDis(K2(i))
            index(i)=K1(i);
        else
            index(i)=K2(i);
        end
    else
        index(i)=K2(i);
    end
    
end
Parents=Population(index,:);
end