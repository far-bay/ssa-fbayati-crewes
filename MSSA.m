clc
close all
clear all
load seismic_color.mat
w = warning('query','last');
id = w.identifier;
warning('off',id);
rmpath('folderthatisnotonpath');
warning('off')
load  Data_3D_2_60; %3D data snr=2 gap=60%
arr="a"; %if a the algorithm uses adaptive rank reduction otherwise uses traditional rank reduction
wrr="w"; %if w the algorithm uses weighted rank reduction otherwise uses traditional rank reduction
orr="0"; %if o the algorithm uses optimal rank reduction otherwise uses traditional rank reduction
method='AWRR'; %depending on the method change to TRR WRR ARR AWRR ORR
figure; wigb(datakill(:,:,1));
Nitr = 20; %number of iteration for the interpolation algorithm to converge the amplitude
dt = 0.002;
dx=10; dy=50;
flag = 1; %if the data is noisy flag=1 if it is noisless flag=0
rank = 12;% select rank depending on the number of linear events

alpha=[];% weighting coeficient for the iterative algorithm
for i=1:Nitr
    alpha(i)=((1-10^(-6))/(1-Nitr))*(i-Nitr);
end
d0 = datakill; d1 = data; d3=addnoise;
%%
[N,n11,n12] = size(data);
f1=1; f2 = 150;
fS = 1/dt; fN = fS/2; fhigh = fS;
if f2>fN; f2=fN; end
N1 = floor(f1*N*dt); N2 = floor(f2*nt*dt);
freqaxis=(0:N);
epsilon = 0.001;
T=datakill ~= 0;
[N,n1,n2] = size(data);
win1 = 23; over1 = floor(win1/2);%spatial window1
win2 = 11; over2 = floor(win2/2);%spatial window2
w1 = win1-over1; nwin1 = floor((n1-over1)/(w1)); mid1 = floor(win1/2+1);
w2 = win2-over2; nwin2 = floor((n2-over2)/(w2)); mid2 = floor(win2/2+1);
mid=mid1*mid2;
%%
index=[];TK = 1;Sob=[];
datakillhed = zeros(N+2,n1,n2);
for i=1:n2
    datakillhed(1,1:n1,i) = 1:n1; %number of receivers for the header
end

for j = 1:n1
    datakillhed(2,j,1:n2) = 1:n2; %number of the shots for the header
end
datakillhed(3:end,1:n1,1:n2)=datakill;
sampling = datakill~=0;
for ii = 1:nwin1
    for jj = 1:nwin2
        datakillw   (:,win1*(ii-1) +1 : win1*ii, win2*(jj-1) +1 : win2*jj ) = datakillhed(:,(ii- 1)*w1+1:ii*w1+over1, (jj- 1)*w2+1:jj*w2+over2) ;
        samplingw   (:,win1*(ii-1) +1 : win1*ii, win2*(jj-1) +1 : win2*jj ) = sampling (:,(ii- 1)*w1+1:ii*w1+over1, (jj- 1)*w2+1:jj*w2+over2);
        %   datakillwOBS(:,win1*(ii-1) +1 : win1*ii, win2*(jj-1) +1 : win2*jj ) = datakillh(:,(ii- 1)*w1+1:ii*w1+over1, (jj- 1)*w2+1:jj*w2+over2);
    end
end
%% function starts
fftdatakill = fft(datakillw(3:end,:,:));%%
figure; wigb(datakillw(3:end,:,1));
DD=datakill ~= 0;
fftdatakill_bp = zeros(size(fftdatakill));
fftdatakill_bp(floor(f1*N*dt)+2:floor(f2*N*dt),:,:) = fftdatakill(floor(f1*N*dt)+2:floor(f2*N*dt),:,:);
tic
for ii = 1:nwin1
    for jj = 1:nwin2
        fftdatakill_recovered   = fftdatakill_bp(:,win1*(ii-1) +1 : win1*ii, win2*(jj-1) +1 : win2*jj );
        fftdatakill0 = fftdatakill_bp(:,win1*(ii-1) +1 : win1*ii, win2*(jj-1) +1 : win2*jj );
        TTT          =      samplingw(:,win1*(ii-1) +1 : win1*ii, win2*(jj-1) +1 : win2*jj );
        %%
        for iter = 1:Nitr
            for f = 1:N
                if f>N1 ||f <N2
                    dd3 = fftdatakill0(f,:,:);
                    TT = TTT(f,:,:);
                    Sob=reshape(dd3,win1,win2);
                    [blockmat1,recovered] = makeblockhankel(fftdatakill_recovered,f,mid1,mid2,win1,win2);
                    blockmat1(isnan(blockmat1))=0;
                    if arr=="a"
                        [U,S1,V,k]=ARR3(blockmat1,rank);
                    else
                        blockmat1(isnan(blockmat1))=0;
                        [U,S,V]=svd(blockmat1);
                        if S(:,:)==0
                            S1=zeros(size(S));
                            k=rank;
                        else
                            k=rank;
                            S1=S;
                            S1(k+1:end,k+1:end)=0;
                        end
                    end
                    if wrr=='w'
                        [U,S1,V]=WRR(blockmat1,mid,k);
                    end
                    if orr=='o'
                        [U,S1,V]=ORR(blockmat1,mid,k);
                    end
                    svdh=U*S1*V';
                    MM(1,1:mid1) = mid2;
                    svdblockcell = mat2cell(svdh,MM,MM);
                    %%
                    for ik = 1:mid1
                        for jk = 1:mid1
                            svdhankvec{ik+jk-1} = svdblockcell{ik,jk};
                        end
                    end
                    for ikk = 1:win1
                        hanksvd = svdhankvec{ikk};
                        % anti diag ave
                        for mmm = 1:mid2
                            if mmm+mid2- 1<= win2; s(mmm) = hanksvd(1,mmm); end
                            if mmm>= 2; s(mmm+mid2- 1) = hanksvd(mmm,mid2); end
                        end
                        h3(ikk,:) = s;
                    end
                    %% Iterative Algorithm
                    h4(:,:) = h3;
                    TT = reshape(TT,win1,win2);
                    iT = ones(size(TT));
                    kkk = iT- TT;
                    %%for clean data
                    if flag==0
                        index1 = norm(((recovered)-(Sob + alpha(iter) * kkk .* h4)),'fro');
                        if index1 <= epsilon; Sob = recovered ; TK=0;
                        else; TK=1; end
                        recovered = Sob + TK * kkk .* h4;
                    end
                    %% for noisy data
                    if flag ==1
                        index1 = norm(((recovered)-(alpha(iter) * Sob + (1-alpha(iter)) * TT .* h4 + kkk .* h4)),'fro');
                        if index1 <= epsilon; Sob = recovered ; TK=0;
                        else; TK=1; end
                        recovered = alpha(iter) * Sob + (1-alpha(iter)) * TT .* h4 + kkk .* h4;
                        %%
                        index(f,iter)=index1;
                        index(f,iter+Nitr)=TK;
                    end
                    h5(f,:,:) = recovered;
                else h5(f,:,:)=0;
                end
            end
            fftdatakill_recovered = h5;
        end
        final(:,win1*(ii-1) +1 : win1*ii, win2*(jj-1) +1 : win2*jj) = fftdatakill_recovered;
    end
end
total_time=toc;
%%
modelhed=zeros(N+2,size(final,2),size(final,3));
model = real(2*(ifft(final)));
modelhed(1:2,:,:) = datakillw(1:2,:,:);
modelhed(3:end,:,:) = model(:,:,:);
[output1]= window_back(modelhed);
%%

output=output1;%(:,floor(aa/2)+1:n11+floor(aa/2),floor(bb/2)+1:n12+floor(bb/2));
a1=1; a2=990;
[N11,nn1,nn2]=size(output);
svdsec1=reshape(output,N11,nn1*nn2);
datakill1=d0(:,1:nn1,1:nn2);
datakill2=reshape(datakill1,N,nn1*nn2);
%%
samplhead=zeros(N+2,size(final,2),size(final,3));
samplhead(1:2,:,:) = datakillw(1:2,:,:);
samplhead(3:end,:,:) = samplingw(:,:,:);
[sample]= window_back(samplhead);
sampelak=[];
samplak=reshape(sample,N11,nn1*nn2);
%%
d0 = d0(:,1:nn1,1:nn2);
d2 = output(3:end,1:nn1,1:nn2);
d1=d1(:,1:nn1,1:nn2);
d3=d3(:,1:nn1,1:nn2);
ISNR = SNR(d0,d1);
OSNR = SNR(d2,d1);
noiseSNR = SNR(d3,d1);
excel(:)=[snr, per, dt, dx, win1, win2, over1, over2, Nitr, arr, wrr, orr, method, total_time, ISNR, OSNR, rank];
%%
figure; imagesc(datakill2(:,a1:a2));set(gca,'FontSize',25);title('Input data');set(gca,'FontSize',25);yticks([0 100 200 300 400 500])
yticklabels({'0' '0.4' '0.8' '1.2' '1.6' '2'});xlabel('Trace Number'); ylabel('Time(s)')
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

data1=d1(:,1:nn1,1:nn2); data2=reshape(data1,N,nn1*nn2);
figure; imagesc(data2(:,a1:a2));set(gca,'FontSize',25);title(' Clean and complete data');yticks([0 100 200 300 400 500])
yticklabels({'0' '0.4' '0.8' '1.2' '1.6' '2'});xlabel('Trace Number'); ylabel('Time(s)')
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

figure;imagesc(svdsec1(3:N11,a1:a2));set(gca,'FontSize',25);title( method );yticks([0 100 200 300 400 500])
yticklabels({'0' '0.4' '0.8' '1.2' '1.6' '2'});xlabel('Trace Number'); ylabel('Time(s)')
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

diff=data2(:,a1:a2)-svdsec1(3:N11,a1:a2);diff(1,:)=max(max(max(datakill)));
figure; imagesc(diff); set(gca,'FontSize',25);title(strcat( method, ' Residuals'));yticks([0 100 200 300 400 500])
yticklabels({'0' '0.4' '0.8' '1.2' '1.6' '2'});xlabel('Trace Number'); ylabel('Time(s)')
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

%%
figure; wigb(datakill2(100:300,400:460));set(gca,'FontSize',15);title('Input data');
xlabel('Trace Number'); ylabel('Time(s)')
yticks([0 100 200 300 400 500]);yticklabels({'0.4' '0.8' '1.2' }),xtick([1 30 60]); xticklabels({'400' '430' '460' })
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

figure; wigb(data2(100:300,400:460));set(gca,'FontSize',15);title(' Clean and complete data');
yticks([0 100 200 300 400 500]);yticklabels({'0.4' '0.8' '1.2' });xlabel('Trace Number'); ylabel('Time(s)')
xticks([1 30 60 ]);  xticklabels({'400' '430' '460'});
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

figure;wigb(svdsec1(100:300,400:460));set(gca,'FontSize',15);title( method );yticks([101 200 300 ])
yticks([0 100 200 300 400 500]);yticklabels({'0.4' '0.8' '1.2' });xlabel('Trace Number'); ylabel('Time(s)')
xticks([1 30 60 ]);  xticklabels({'400' '430' '460'});
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

diff=data2(:,a1:a2)-svdsec1(3:N11,a1:a2);diff1=diff(100:300,400:460);diff1(1,:)=max(max(max(datakill)));
figure; wigb(diff1); set(gca,'FontSize',15);title(strcat( method, ' Residuals'));xlabel('Trace Number'); ylabel('Time(s)')
yticklabels({'0.4' '0.8' '1.2'});
yticks([0 100 200 300 400 500]);xticks([1 30 60 ]);  xticklabels({'400' '430' '460'});
xlabel('Trace Number'); ylabel('Time(s)')
colormap(seismic_color);caxis([-1-1/snr 1+1/snr])

