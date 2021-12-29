function [U,S1,V,k]=ARR2(Hankel,rank)
[U, S, V]=svd(Hankel);
if S(:,:)==0
    S1=zeros(size(S));
    k=rank;
else
    if S(2,2)==0; k=rank;
    else
        A=((S(1:end-1,1:end-1).^2./S(2:end,2:end).^2));
        B=A(1:end-3,1:end-3);
        c=diag(B);
        BA=(max(c));
        k=find(c==BA)+1;
       
    end
    if k> size(S,2)
        k=rank;
    end
    S1=S;
    S1(k+1:end,k+1:end)=0;
end
end