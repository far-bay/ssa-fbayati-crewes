
function [blockmat1,recovered] = makeblockhankel(fftdatakill_recovered,f,mid1,mid2,lent1,lent2)
for len = 1:lent1
    recovered = fftdatakill_recovered(f,:,:);
    recovered = recovered(:,1:lent1,1:lent2);
    recovered = reshape(recovered,lent1,lent2);
    d = fftdatakill_recovered(f,len,:);
    h = hankel(d);
    I = h(1:mid2,1:mid2);
    kk{len} = I;
end
for i = 1:mid1
    for j = 1:mid1
        cellhankel{i,j} = kk{i+j-1};
    end
end
blockmat1 = cell2mat(cellhankel);
end
