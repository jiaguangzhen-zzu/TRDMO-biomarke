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
