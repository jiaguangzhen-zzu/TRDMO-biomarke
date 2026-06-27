function [outputArg1,outputArg2] = MMPDENB_PSL(inputArg1,inputArg2)
%% MOEA-PSL identify edge-network biomarkers
%   此处显示详细说明
data=load('E:\数据\TCGA1\Data11_TCGA-KICH.mat');

    experiment_num=30;
    Problem.maxFE=30000;
    Problem.N=300;
    Problem.lower=0;
    Problem.upper=1;
    TPOP = [];
    for fnum=1:4 %length(Sample_net)%邻接矩阵
        fieldName = ['N', num2str(fnum)];
        test_adjacency = data.(fieldName);
       
        Problem.D=size(test_adjacency,1); %维度
        test_net=zeros(Problem.D,Problem.D);
        test_net(test_adjacency~=0)=1;
        %% paired edge :g
        Cons=CONS(test_net);
        Cons(all(Cons==0,2),:)=[];%删除全零行
%         cons=Cons;%约束条件
        Cnum1=size(Cons,1);
        %Cnum2=size(Cons,2);
        g=cell(Cnum1,1);
        %validIndices = true(Cnum2, 1);  % 逻辑索引，用于指示哪些行是维度一致的
        for i=1:Cnum1
            g{i}=find(Cons(i,:));
            %if i > 1 && numel(g{i}) ~= numel(g{i - 1})
                %validIndices(i) = false;  % 标记维度不一致的行
            %end
        end
        %g = g(validIndices, :);  % 保留维度一致的行
        Cnum=size(g,1);
        clear Cons
        %% 计算每个维度的目标函数值
        Dimension=eye(Problem.D,Problem.D);
        CV_num=Calcons(Problem.D,Cnum,g,Dimension);
        D_score=abs(CV_num-Cnum);
        Non_dominated_sol=cell(30,1);
        R=cell(30,1);
        for EXP_NUM=1:experiment_num
            CalNum=0;
            
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
            %             Mask = UniformPoint(Problem.N,Problem.D,'Latin') > 0.5;
             Population = Dec.*Mask;
             [Population,calnum]  = SOLUTION(Population,test_adjacency,g,D_score);
             CalNum=CalNum+calnum;
            
             [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);%%里面的非支配排序要换成自己的
           else
              Dec = ones(Problem.N,Problem.D); 
              
                  numPopRows = size(TPOP, 1);  % 获取 pop 矩阵的行数
            
                 maxPopRows = min(100, size(TPOP, 1));  % 从 pop 中取出的最大行数为100

                  if size(pop, 1) < 100
                     
                         neededRows = 100 - size(TPOP, 1);
                        part1 = TPOP;  
                        part2 = creatpop(neededRows,Problem.D,D_score,g); % 生成补足的随机数据
                  else
                        
                       part1 = TPOP(1:100, :);
                        part2 = [];  % 不需要额外的补足行
                  end 
                  part3 = creatpop(200,Problem.D,D_score,g);

                 Mask = [part1; part2; part3];  % 总共300行
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
                    %MatingPool = TournamentSelection_hamming(Problem.N,Problem.N,Popmarks,FrontNo,CrowdDis);
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
                    MatingPool = TournamentSelection_hamming(Problem.N,Problem.N,Popmarks,FrontNo,CrowdDis);
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
            outputpop=Pop(FrontNo==1,:);
            Non_dominated_sol{EXP_NUM}=outputpop;
            R{EXP_NUM}=RHO;
        end
        
        Non_dominated_sol1=cell2mat(Non_dominated_sol);
        
        [pop,~,~]=unique(Non_dominated_sol1,'rows');%将得到的非支配解进行去重
        
        [functionvalue,~] = Calfunctionvalue(pop,test_adjacency);%计算去重后的非支配解的目标函数
              
            %folderPath = 'E:\迁移种群\目标值'; % 修改为您的文件夹路径
            
            %filename = fullfile(folderPath, 'matrixA.mat');

         
            %save(filename, 'functionvalue');

        
        [FrontNo,~] = NDSort(functionvalue,size(functionvalue,1));%对这些解进行pareto等级排序
        
        POP=pop(FrontNo==1,:);%输出最终的非支配解对应个体
        
        FV=functionvalue(FrontNo==1,:);
        
        [FVV,mod_position]=sortrows(FV);
        
        Final_Pop=POP(mod_position,:);% 输出最终排序过的个体
     
       
        Tpop = vertcat(TPOP, pop);
        [Tpop,~,~]=unique(Non_dominated_sol1,'rows');
        [functionvalue1,~] = Calfunctionvalue(Tpop,test_adjacency);
        [FrontNo1,~] = NDSort(functionvalue1,size(functionvalue1,1));
         [~, idx] = sort(FrontNo1);  % 获取排序后的索引
        tpop = Tpop(idx, :);  % 对 Pop 按照 FrontNo 的值进行排序

     % 取排序后的前 100 个个体作为输出种群
      if size(tpop, 1) < 100
          TPOP = tpop ;
      else TPOP = tpop(1:100, :);
      end

        %% 统计多模态出现的次数
        FVV(:,2)=-FVV(:,2);
        A=tabulate(FVV(:,1));
        B= A(:,2)~=0;
        c=A(B,:);
        %% 储存数据
            
        
        filename=strcat('E:\结果\Data10\TPRORBM\TPRORBM1',num2str(fnum),'_MMPDENB-RBM_NDSol.mat');
        save(filename,'FVV');
        filename1=strcat('E:\结果\Data10\TPRORBM\TPRORBM1',num2str(fnum),'_MMPDENB-RBM_Num.mat');
        save(filename1,'c');
        filename2=strcat('E:\结果\Data10\TPRORBM\TPRORBM1',num2str(fnum),'_MMPDENB-RBM_POP.mat');
        save(filename2,'Final_Pop');
        filename3=strcat('E:\结果\Data10\TPRORBM\TPRORBM1',num2str(fnum),'_MMPDENB-RBM_boxchart.mat');
        save(filename3,'Non_dominated_sol');
        filename4=strcat('E:\结果\Data10\TPRORBM\TPRORBM1',num2str(fnum),'_MMPDENB-RBM_pho.mat');
        save(filename4,'R');

    end
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