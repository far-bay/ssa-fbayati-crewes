function [output]= window_back(inputdata)

for i=1: size(inputdata,3)
    BB= inputdata(:,:,i);
    CC=reshape(BB,size(BB,1),size(BB,2));
    uniqueNEW1 = unique(CC(1,:));
    uniqNEW(1,:,i) = uniqueNEW1;
end
CC=[];
BB=[];
for i=1: size(inputdata,2)
    BB= inputdata(:,i,:);
    CC=reshape(BB,size(BB,1),size(BB,3));
    uniqueinputdata = unique(CC(2,:));
    uniqNEW3(1,i,:) = uniqueinputdata;
end
[N1,nn]=size(uniqueNEW1);
[N1,nn2]=size(uniqueinputdata);
for j=1:size(inputdata,3)
    for i=1: nn
        AA=inputdata(:,:,j);
        AA=reshape(AA,size(inputdata,1),size(inputdata,2));
        pp = find((AA(1,:))== uniqNEW(i));
        lp = length(AA(pp));
        sumA = sum(AA(:,pp),2);
        normstck = sumA/lp;
        final(:,i,j) = normstck;
    end
end
for j=1:size(final,2)
    for i=1: nn2
        BB=final(:,j,:);
        BB=reshape(BB,size(inputdata,1),size(inputdata,3));
        uniq=uniqNEW3(1,i,:);
        pp = find((BB(2,:))== uniq(i));
        lp = length(BB(pp));
        sumA = sum(BB(:,pp),2);
        normstck = sumA/lp;
        final1(:,j,i) = normstck;
    end
end
output=final1;
 end