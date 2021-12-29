function [U,S1,V,S,Sr]=WRR(Hankel,mid,rank)
[U, S, V]=svd(Hankel);
if S(:,:)==0;
    S1=zeros(size(S));
else
    I=eye(mid-rank);
    ss=S(rank+1:end,rank+1:end);
    for i=1:rank
        term=S(i,i)^2.*I-ss.^2;
        d=((1/mid)*(trace(S(i,i)*(term)^(-1)))).^2;
        dr=(2/mid^2)*(trace(S(i,i)*(term)^(-1)))*(trace((term).^(-1)-2*S(i,i)^2*(term)^(-2)));
        W(i,i)=(-2/S(i,i).*(d./dr));
    end
    WR=zeros(mid);
    Sr=zeros(mid);
    Sr(1:rank,1:rank)=S(1:rank,1:rank);
    WR(1:rank,1:rank)=W(1:rank,1:rank);
    S1=WR*Sr;
    
end
% result=U*WR*S*V';

end
