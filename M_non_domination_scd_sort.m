function [FrontNo,SpCrowdDis] = M_non_domination_scd_sort(Population,functionvalue)

 %% 先除去相同的个体
 

     [functionvalue_sort,ia,~] = unique(functionvalue,'rows');  %% C = A(ia)并且A = C(ic)

    
    fnum=0;
    V=2;
     frontvalue=zeros(1,size(functionvalue,1));
    cz=zeros(1,size(functionvalue_sort,1));
    s=zeros(size(functionvalue_sort,1),1);
    %% 给每个个体分配等级
    while ~all(cz)
        fnum=fnum+1;
        d=cz;
        for i=1:size(functionvalue_sort,1)
            s(i,:)=0;
            if ~d(i)
                for j=i+1:size(functionvalue_sort,1)
                    if ~d(j)
                        for k=2:V
                            if functionvalue_sort(i,k)<=functionvalue_sort(j,k)
                                d(j)=1;
                            else
                                d(j)=0;
                            end
                        end
                    end
                end
                frontvalue(:,ia(i))=fnum;
                cz(i)=1;
            end
        end
    end

    %% 给多模态解分配前沿编号   
    dis=pdist2(functionvalue,functionvalue);
    while ~all(frontvalue)              %% 如果前沿面有0则循环

          Multimodal=find(frontvalue==0);  %% 找到多模态解的位置
          Multimodal_position=find(dis(Multimodal(1),:)==0);  %%找到目标函数中与第一个多模态解目标函数一样解的位置
          front=frontvalue(Multimodal_position);          %% 找出多模态解的前沿面
          have_fornt_position=find(front~=0);               %% 找出已经赋值前沿面的位置
          Multimodal_position(have_fornt_position)=[];        %% 删除该位置
          frontvalue(Multimodal_position)=front(have_fornt_position);  %% 对没分配前沿面的分配前沿面

    end
    
    FrontNo= frontvalue;
    
    %% Calculate the special crowding distance of each solution
    SpCrowdDis_Obj = ModifiedCrowdingDistance(functionvalue,FrontNo);
     SpCrowdDis_Dec = ModifiedCrowdingDistance_dec(Population,FrontNo,functionvalue);
    %% Limit the size of Population 
     SpCrowdDis = max(SpCrowdDis_Obj,SpCrowdDis_Dec);
%     SpCrowdDis = SpCrowdDis_Obj;
end

function CrowdDis = ModifiedCrowdingDistance(PopObj,FrontNo)
    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Sol_in_Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Sol_in_Front,:),[],1);
        Fmin  = min(PopObj(Sol_in_Front,:),[],1);
        for i = 1 : M    %%拥挤距离是两个目标函数的拥挤距离累加和
            [~,Rank] = sortrows(PopObj(Sol_in_Front,i));
            CrowdDis(Sol_in_Front(Rank(1))) = CrowdDis(Sol_in_Front(Rank(1))) + 1;
            for j = 2 : length(Sol_in_Front)-1
                if Fmax(i) == Fmin(i)
                    CrowdDis(Sol_in_Front(Rank(j))) = CrowdDis(Sol_in_Front(Rank(j)))+1;
                else
                    CrowdDis(Sol_in_Front(Rank(j))) = CrowdDis(Sol_in_Front(Rank(j)))+(PopObj(Sol_in_Front(Rank(j+1)),i)-PopObj(Sol_in_Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
                end
            end
        end
    end
end
function CrowdDis = ModifiedCrowdingDistance_dec(Pop_dec,FrontNo,functionvalue)
Pop_dec=double(Pop_dec);
    N    = size(Pop_dec,1);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Sol_in_Front = find(FrontNo==Fronts(f));

            [~,Rank] = sortrows(functionvalue(Sol_in_Front,1));   %% 对前沿面上的解按照目标函数1从小到大排序
            if length(Sol_in_Front)-1==0
                 CrowdDis(Sol_in_Front(Rank(1))) = CrowdDis(Sol_in_Front(Rank(1)))+1;
            else
                CrowdDis(Sol_in_Front(Rank(1))) = CrowdDis(Sol_in_Front(Rank(1))) + 1;
                CrowdDis(Sol_in_Front(Rank(end))) = CrowdDis(Sol_in_Front(Rank(end))) + 1;
                for j = 2 : length(Sol_in_Front)-1
                    CrowdDis(Sol_in_Front(Rank(j))) = (pdist2(Pop_dec(Sol_in_Front(Rank(j)),:),Pop_dec(Sol_in_Front(Rank(j-1)),:),'hamming')+ ...
                                                         pdist2(Pop_dec(Sol_in_Front(Rank(j)),:),Pop_dec(Sol_in_Front(Rank(j+1)),:),'hamming'))/2;
                end
            end
    end
end
