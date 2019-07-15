### Datapreprocess

```
clear
load feature.mat EntrezGene RNAtype HugoSymbol Length
load labels.mat sampleID CancerType
load expdata.mat exp sampleid features

%label
[samples,isID,isid]=intersect(sampleID,sampleid);
labels=CancerType(isID);
po1=find((ismember(labels,{'HCC'})~=0)~=0);
po2=find((ismember(labels,{'CRC'})~=0)~=0);
po3=find((ismember(labels,{'PAAD'})~=0)~=0);
ne=find((ismember(labels,'Healthy')~=0)~=0);
%% 训练集
m=0.7;%正样本百分数
n=0.8;%负样本百分数
nh1=ceil(length(po1)*m);%hcc正个数
nh2=ceil(length(po2)*m);%crc正个数
nh3=ceil(length(po3)*m);%paad正个数
nn1=ceil(length(ne)*n);%正常样本数
c11=po1(randperm(numel(po1)));%随机取
c12=po2(randperm(numel(po2)));%随机取
c13=po3(randperm(numel(po3)));%随机取
c22=ne(randperm(numel(ne)));%随机取
ch1=c11(1:nh1);
ch2=c12(1:nh2);
ch3=c13(1:nh3);
cn=c22(1:nn1);
trpo=[ch1;ch2;ch3];
trne=cn;
train=[exp(ch1,:);exp(ch2,:);exp(ch3,:);exp(cn,:)];
trainid=[sampleid(ch1,:);sampleid(ch2,:);sampleid(ch3,:);sampleid(cn,:)];

%% 测试集
th1=nh1+1;%测试集正样本数起始计数点
th2=nh2+1;
th3=nh3+1;
tn2=nn1+1;%测试集负样本数起始计数点
tc1=c11(th1:end);
tc2=c12(th2:end);
tc3=c13(th3:end);
nc1=c22(tn2:end);
tepo=[tc1;tc2;tc3];
tene=nc1;
test=[exp(tc1,:);exp(tc2,:);exp(tc3,:);exp(nc1,:)];
testid=[sampleid(tc1,:);sampleid(tc2,:);sampleid(tc3,:);sampleid(nc1,:)];
%% 验证集

%% Mann-Whitney U test ：
pm=[];hm=[];
[row,col]=size(train);
trainpo=exp(trpo,:);trainne=exp(trne,:);
testpo=exp(tepo,:);testne=exp(tene,:);
for u=1:col;
    [pm2,hm2]=ranksum(trainpo(:,u),trainne(:,u));
    pm=[pm pm2];
    hm=[hm hm2];
end
isone=find(hm==1);
trainTs=trainpo(:,isone);
trainNs=trainne(:,isone);
trainGs=features(isone,1);
trainX=[trainTs;trainNs];
trainY=[ones(length(trpo),1);-ones(length(trne),1)];
trainID=trainid;
save Datatotrain.mat trainX trainY trainGs trainID
testTs=testpo(:,isone);
testNs=testne(:,isone);
testGs=features(isone,1);
testX=[testTs;testNs];
testY=[ones(length(tepo),1);-ones(length(tene),1)];
testID=testid;
save Datatotest.mat test testX testY testGs testID
```

### coeprediction

```
clear
load trycoef.mat coe fea allacc
load Datatotrain.mat trainX trainY trainID trainGs
load Datatotest.mat testX testY testID testGs

unigene=trainGs;
ACCG=[];
ACCT=[];
xcal2=trainX;
tscore2=testX;

for ind=1:200;
secoef=coe(:,ind);
sefea=fea(ind,:);
sefea(202-ind:end)=[];
secoef(203-ind:end)=[];
segene=unigene(sefea);

[lia,loc]=ismember(segene,unigene);
xcal2=xcal2(:,loc);
tcal=testX(:,loc);
[xcal2,para1,para2]=pretreat(xcal2,'autoscaling');
tcal=pretreat(tcal,'autoscaling',para1,para2);

%% predict
gpre=xcal2*secoef(1:end-1)+secoef(end);
yg=sign(gpre);
tpre=tcal*secoef(1:end-1)+secoef(end);
yt=sign(tpre);

%train
nh1=length(find(trainY==1)==1);
nn1=nh1+1;
yg1=yg(1:nh1);
yg1s=length(find(yg1==1));%sensitive TP
yg1r=length(find(yg1==-1));%resistant FN
yg2=yg(nn1:end);
yg2s=length(find(yg2==1));%sensitive FP
yg2r=length(find(yg2==-1));%resistant TN

sn=yg1s/(yg1s+yg1r);
sp=yg2r/(yg2r+yg2s);
acc=(yg1s+yg2r)/(yg1s+yg1r+yg2s+yg2r);
mcc=(yg1s*yg2r-yg1r*yg2r)/sqrt((yg1s+yg2s)*(yg1s+yg1r)*(yg2r+yg2s)*(yg1r+yg2r));
ytrain=[sn;sp;acc;mcc];
%test
nh2=length(find(testY==1)==1);
nn2=nh2+1;
yt1=yt(1:nh2);
yt1s=length(find(yt1==1));%sensitive
yt1r=length(find(yt1==-1));%resistant
yt2=yt(nn2:end);
yt2s=length(find(yt2==1));%sensitive
yt2r=length(find(yt2==-1));%resistant

sn=yt1s/(yt1s+yt1r);
sp=yt2r/(yt2r+yt2s);
acc=(yt1s+yt2r)/(yt1s+yt1r+yt2s+yt2r);
mcc=(yt1s*yt2r-yt1r*yt2r)/sqrt((yt1s+yt2s)*(yt1s+yt1r)*(yt2r+yt2s)*(yt1r+yt2r));
ytest=[sn;sp;acc;mcc];

ACCG=[ACCG ytrain];
ACCT=[ACCT ytest];


xcal2=trainX;
tscore2=testX;
para1=[];
para2=[];
end
plotACC=[ACCG ACCT];
ACC=[ACCG;ACCT];
for i=1:2
figure (2)
subplot(2,1,i)
for j=1:4
plot(plotACC(j,(200*(i-1)+1:200*i)))
set(gca,'xtick',[0:10:200],'xticklabel',[200:-10:0])
hold on 
end
box off
legend('sn','sp','acc','mcc','location','southwest')
end

n=178; 
sefea=fea(n,:);
sefea(202-n:end)=[];
secoef(203-n:end)=[];
segene=unigene(sefea);
accse=[ACCG(:,n) ACCT(:,n)];
save segene62.mat segene secoef accse;
```

