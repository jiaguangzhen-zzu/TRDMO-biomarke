function PopObj        = extract_d(Population)
%% PopObj=Population.decs
        PopObj=zeros(size(Population,2),size(Population{1, 1}.dec,2));
        for i=1:size(PopObj,1)
            PopObj(i,:)=Population{1,i}.dec;
        end
end