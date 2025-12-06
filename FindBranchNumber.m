function Output=FindBranchNumber(Input,branch)
row = size(Input,1);
for i = 1:row
    if Input(i,1) > Input(i,2)
        a = Input(i,1);
        Input(i,1) = Input(i,2);
        Input(i,2) = a;
        m(i) = -1;
    else
        m(i) = 1;
    end   
end

for i = 1:row
    for j = 1:size(branch,1)
        if Input(i,1) == branch(j,1)
            if Input(i,2) == branch(j,2)
                num(i) = j;
            end
        end
    end
end

Output = [m;num]';
end