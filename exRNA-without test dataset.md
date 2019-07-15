## exRNA

### classification without test dataset

#### get data
```
fileID = fopen('feature_name.txt','r');
slCharacterEncoding('ISO-8859-1')
nNumberCols =4;
format =  repmat('%s', [1 nNumberCols]);
A=textscan(fileID,format,109441,'Headerlines',1,'Delimiter','\t');
EntrezGene=A{1,1};
RNAtype=A{1,2};
HugoSymbol=A{1,3};
Length=A{1,4};
fclose(fileID);
save feature.mat EntrezGene RNAtype HugoSymbol Length
%% get labels
clear
fileID = fopen('sample_classes.txt','r');
slCharacterEncoding('ISO-8859-1')
nNumberCols =2;
format =  repmat('%s', [1 nNumberCols]);
A=textscan(fileID,format,86,'Headerlines',1,'Delimiter','\t');
sampleID=A{1,1};
CancerType=A{1,2};
fclose(fileID);
save labels.mat sampleID CancerType
%% get expression
clear
[exp,text]=xlsread('featurecounts.xlsx','featurecounts');
[row,col]=size(text);
features=text((2:row),1);
sampleid=text(1,(2:col));
sampleid=sampleid';
exp=exp';
save expdata.mat exp sampleid features
```

#### datapreprocess

```
clear
load feature.mat EntrezGene RNAtype HugoSymbol Length
load labels.mat sampleID CancerType
load expdata.mat exp sampleid features

%label
[samples,isID,isid]=intersect(sampleID,sampleid);
labels=CancerType(isID);
positive=find((ismember(labels,{'HCC','CRC','PAAD'})~=0)~=0);
negative=find((ismember(labels,'Healthy')~=0)~=0);
expT=exp(positive,:);
expN=exp(negative,:);
labelT=ones(length(positive),1);
labelN=-ones(length(negative),1);
sampleT=sampleid(positive);
sampleN=sampleid(negative);
%keyboard
%% Mann-Whitney U test ：
pm=[];hm=[];
for u=1:length(HugoSymbol)
    [pm2,hm2]=ranksum(expT(:,u),expN(:,u));
    pm=[pm pm2];
    hm=[hm hm2];
end
isone=find(hm==1);
trainTs=expT(:,isone);
trainNs=expN(:,isone);
trainGs=HugoSymbol(isone);
features=features(isone,1);
trainX=[trainTs;trainNs];
trainY=[labelT;labelN];
trainID=[sampleT;sampleN];
save Datatotrain.mat trainX trainY trainGs features trainID
```

#### PLS training

```
clear
load Datatotrain.mat trainX trainY trainGs features trainID

geneid=1:length(trainGs);
geneid2=geneid;
Xcal=trainX;
ycal=trainY;
sn=[];
sp=[];
acc=[];
auc=[];
mcc=[];
coef2=[];
B2=[];
coe=[];
coef3=[];
VIP2=[];
xh=10;%循环次数

d1=mod(length(trainGs),100)+100;%mod取余数
for j=1:xh
A=5;%The maximal number of LVs for cross-validation
K=3;%fold. when K = m, it is leave-one-out CV
LDA=plslda(Xcal,ycal,A,K,'autoscaling',0,0);
B2=[B2 LDA.B];%重要性指标
VIP2=[VIP2 LDA.VIP];
coef3=[coef3 LDA.coef];%权重系数
end
coef2=median(coef3,2);
B=median(B2,2);
[~,index]=sort(abs(B));
VIP=median(VIP2,2);
[~,index]=sort(VIP);
index(1:d1)=[];
Xcal=Xcal(:,index);
geneid=geneid(:,index);
B2=[];
VIP2=[];
coef3=[];
coef2=[];
clear B index VIP

fea=geneid;
d2=fix((length(geneid)-200)/200);%fix:截尾取整
for i=1:d2
for j=1:xh
A=5;%The maximal number of LVs for cross-validation
K=3;%fold. when K = m, it is leave-one-out CV
LDA=plslda(Xcal,ycal,A,K,'autoscaling',0,0);
B2=[B2 LDA.B];%重要性指标
VIP2=[VIP2 LDA.VIP];
coef3=[coef3 LDA.coef];%权重系数
end
coef2=median(coef3,2);
B=median(B2,2);
[~,index]=sort(abs(B));
VIP=median(VIP2,2);
[~,index]=sort(VIP);
index(1:200)=[];
Xcal=Xcal(:,index);
geneid=geneid(:,index);
B2=[];
VIP2=[];
coef3=[];
coef2=[];
end
clear B index VIP

fea=geneid;
for i=1:200
for j=1:xh
A=5;%The maximal number of LVs for cross-validation
K=3;%fold. when K = m, it is leave-one-out CV
LDA=plslda(Xcal,ycal,A,K,'autoscaling',0,0);
B2=[B2 LDA.B];%重要性指标
VIP2=[VIP2 LDA.VIP];
coef3=[coef3 LDA.coef];%权重系数
end
coef2=median(coef3,2);
B=median(B2,2);
[~,index]=sort(abs(B));
VIP=median(VIP2,2);
[~,index]=sort(VIP);
index(1:1)=[];
Xcal=Xcal(:,index);
geneid=geneid(:,index);
feature=[geneid zeros(1,i)];
fea=[fea;feature];%[100个特征：99+0个特征]
se2=LDA.sensitivity;
sp2=LDA.specificity;
acc2=1-LDA.error;
mcc2=LDA.MCC;
auc2=LDA.AUC;
coef2=[coef2;zeros((i-1),1)];%[101个系数；0个系数]
coe=[coe coef2];
%keyboard
sn=[sn se2];
sp=[sp sp2];
%ave2=[sn;sp];
acc=[acc acc2];
mcc=[mcc mcc2];
auc=[auc auc2];
%ave=mean(ave2);
B2=[];
VIP2=[];
coef3=[];
end
allacc=[sn;sp;acc;mcc];
save trycoef.mat fea coe allacc

figure (1)
for i=1:4
    plot (allacc(i,:))
    hold on
end
box off
```

#### coeprediction

```
clear
load trycoef.mat coe fea allacc
load Datatotrain.mat trainX trainY trainGs features trainID

unigene=trainGs;
ACCG=[];
xcal2=trainX;

for ind=181;
secoef=coe(:,ind);
sefea=fea(ind,:);
sefea(202-ind:end)=[];
secoef(203-ind:end)=[];
segene=unigene(sefea);

[lia,loc]=ismember(segene,unigene);
xcal2=xcal2(:,loc);
[xcal2,para1,para2]=pretreat(xcal2,'autoscaling');
keyboard
%% predict
gpre=xcal2*secoef(1:end-1)+secoef(end);
yg=sign(gpre);

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
ACCG=[ACCG ytrain];

xcal2=trainX;
para1=[];
para2=[];
end
plotACC=[ACCG];
ACC=[ACCG];
figure (2)
for j=1:4
plot(plotACC(j,(1:200)))
set(gca,'xtick',[0:10:200],'xticklabel',[200:-10:0])
hold on 
end
box off
legend('sn','sp','acc','mcc','location','southwest')

n=181; 
sefea=fea(n,:);
sefea(202-n:end)=[];
secoef=coe(:,n);
secoef(203-n:end)=[];
segene=unigene(sefea);
accse=ACCG(:,n);
save segene12.mat segene secoef accse ;
```

------

### function scripts

#### vip

```
function VIP=vip(X,y,T,W,Q)
%+++ Calculate the vip for each variable to the response;
%+++ T,W,Q are output from plsnipals.m
%+++ T: socres, which can be obtained by pls_nipals.m 
%+++ W: weight, which can be obtained by pls_nipals.m 
%+++ Q: Y-loadings
%+++ VIP=sqrt(p*q/s);
%+++      where, p is a constant, denoting the number of variables in x
%                q stands for the explained variance of Y by each variable
%                s represents the total variance explained by A components
%+++ Reference: Tahir Mehmood et al, Chemometrics and Intelligent Laboratory Systems 118 (2012)62?69
%+++ HDLi


s=diag(T'*T*Q*Q');
%initializing
[m,p]=size(X);
[m,h]=size(T);
%+++ calculate VIP;
VIP=[];
for i=1:p
    weight=[];
    for j=1:h
        weight(j,1)= (W(i,j)/norm(W(:,j)))^2;
    end
    q=s'*weight;  % explained variance by variable i
    VIP(i)=sqrt(p*q/sum(s));   
end

%+++ 
```

#### roccurve

```
function F=roccurve(ypred,yreal,flag)
%+++ yreal: with elements 1 or -1;
%+++ ypred: real values.
%+++ flag: 1: plot
%          0: no plot.
%+++ June 11,2008.

if nargin<3;flag=1;end;

yreal=sign(yreal);
thitamin=min(ypred);thitamax=max(ypred);
K=128;
if K>size(yreal,1);K=size(yreal,1);end

thita=linspace(thitamin-0.000001,thitamax+0.000001,K);
Result=zeros(K,2);
for i=1:K
  r=sesp(ypred-thita(i),yreal);
  Result(i,:)=[r.specificity r.sensitivity];    
end
auc=abs(trapz(Result(:,1),Result(:,2)));
if flag==1
  plot(1-Result(:,1),Result(:,2));
  xlabel('1-specificity');ylabel('sensitivity');
end
r=sesp(ypred,yreal);
%+++ OUTPUT
F.value=Result;
F.MCC=r.mcc;
F.sensitivity=r.sensitivity;
F.specificity=r.specificity;
F.accuracy=r.accuracy;
F.AUC=auc; % area under curve.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% subfunction   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=sesp(ypred_value,yreal_binary)
%+++ To calculate the sensitivity(Se) and the specificity(Sp) of 
%+++ a binary classification problem.
%+++ yreal_binary has to be 1 or -1.
%+++ Hongdong Li, Apr.29,2008.

ypred_value=sign(ypred_value);
yreal_binary=sign(yreal_binary);
p=0;n=0;o=0;u=0;LEN=length(yreal_binary);
for i=1:LEN
    if  yreal_binary(i)==1
       if ypred_value(i)==1;p=p+1;else;o=o+1;end
    elseif yreal_binary(i)~=1
       if ypred_value(i)~=1; n=n+1;else;u=u+1;end
    end   
end

%+++ output
result.sensitivity=p/(p+o);
result.specificity=n/(n+u);
result.mcc=(p*n-u*o)/sqrt((p+u)*(p+o)*(n+u)*(n+o));
result.accuracy=(p+n)/LEN;
```

#### rescaling

```
function [X,para1,para2]=rescaling(X,method,para1,para2)
%+++   data pretreatment
%+++ HD Li, Central South University


if nargin==2
  [Mx,Nx]=size(X);
   if strcmp(method,'autoscaling')
    para1=mean(X);para2=std(X);
   elseif strcmp(method,'center')
    para1=mean(X);para2=ones(1,Nx);
   elseif strcmp(method,'unilength')
    para1=mean(X);
    for j=1:size(X,2);
    para2(1,j)=norm(X(:,j)-para1(j));
    end
   elseif strcmp(method,'minmax')
    para1=min(X);maxv=max(X);
    para2=maxv-para1;  
   elseif strcmp(method,'pareto');
    para1=mean(X);para2=sqrt(std(X));
   else
    display('Wrong data pretreat method!');
   end
   
   for i=1:Nx
      if para2(i)==0
          X(:,i)=zeros(Mx,1);
      else
     X(:,i)=(X(:,i)-para1(i))/para2(i);
      end
   end
   
elseif nargin==4
   [Mx,Nx]=size(X);
   for i=1:Nx    
       if para2(i)==0
          X(:,i)=zeros(Mx,1);
      else
     X(:,i)=(X(:,i)-para1(i))/para2(i);
       end
   end
end
```

#### repeatEntries

```
function out = repeatEntries(val,kTimes)
%REPEATENTRIES fills a matrix with k repeats the rows of the input matrix
%
% SYNOPSIS out = repeatEntries(val,kTimes)
%
% INPUT    val    : matrix (or vectors) containing the rows to repeat (works for strings, too)
%          kTimes : number of repeats of each row (scalar or vector of size(vlaues,1))
%
% OUTPUT   out    : matrix of size [sum(kTimes) size(values,2)] containing
%                   repeated entries specified with k
%
% EXAMPLES     repeatEntries([1;2;3;4],[2;3;1;1]) returns [1;1;2;2;2;3;4]
%
%              repeatEntries([1;2;3;4],2) returns [1;1;2;2;3;3;4;4]
%
% c: jonas, 2/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===========
% test input
%===========

% nargin
if nargin ~= 2 || isempty(val) || isempty(kTimes)
    error('two non-empty input arguments are needed!')
end

% size
valSize = size(val);
if length(valSize)>2
    error('only 2D arrays supported for val')
end



% decide whether we have scalar k
numK = length(kTimes);
if numK == 1
    scalarK = 1;
elseif numK ~= valSize(1)
    error('vector k must have the same length as the number of rows in val or be a scalar')
else
    % check again whether we could use scalar k
    if all(kTimes(1) == kTimes)
        scalarK = 1;
        kTimes = kTimes(1);
    else
        scalarK = 0;
    end
end

% do not care about size of k: we want to make a col vector out of it - and
% this vector should only contain nonzero positive integers
kTimes = round(kTimes(:));
% if there are any negative values or zeros, remove the entry
if scalarK && kTimes < 1
    out = [];
    return
end
if ~scalarK
    badK = kTimes < 1;
    kTimes(badK) = [];
    val(badK,:) = [];
    % update valSize
    valSize = size(val);
    if any(valSize==0)
        out = [];
        return
    end
end
%kTimes = max(kTimes,ones(size(kTimes)));


%============
% fill in out
%============

% first the elegant case: scalar k
if scalarK

    % build repeat index matrix idxMat
    idxMat = meshgrid( 1:valSize(1), 1:kTimes(1) );
    idxMat = idxMat(:); % returns [1;1...2;2;... etc]

    out = val(idxMat,:);

    % second: the loop
else

    % init out, init counter
    if iscell(val)
        out = cell(sum(kTimes) , valSize(2));
    else
    out = zeros( sum(kTimes), valSize(2) );
    end
    endct = 0;

    if valSize(2) == 1

        % vector: fill directly

        % loop and fill
        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            out(startct:endct,:) = val(i);
        end % for i=1:valSize(1)

    else

        % matrix: fill via index list

        idxMat = zeros(sum(kTimes),1);

        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            idxMat(startct:endct) = i;
        end % for i=1:valSize(1)
        out = val(idxMat,:);

    end

    % check for strings and transform if necessary
    if ischar(val)
        out = char(out);
    end

end % if doScalar
```

#### pretreat

```
function [X,para1,para2]=pretreat(X,method,para1,para2)
%+++   data pretreatment
%+++ HD Li, Central South University


if nargin==2
  [Mx,Nx]=size(X);
   if strcmp(method,'autoscaling')
    para1=mean(X);para2=std(X);
   elseif strcmp(method,'center')
    para1=mean(X);para2=ones(1,Nx);
   elseif strcmp(method,'unilength')
    para1=mean(X);
    for j=1:size(X,2);
    para2(1,j)=norm(X(:,j)-para1(j));
    end
   elseif strcmp(method,'minmax')
    para1=min(X);maxv=max(X);
    para2=maxv-para1;  
   elseif strcmp(method,'pareto');
    para1=mean(X);para2=sqrt(std(X));
   else
    display('Wrong data pretreat method!');
   end
   
   for i=1:Nx
     X(:,i)=(X(:,i)-para1(i))/para2(i);
   end
   
elseif nargin==4
   [Mx,Nx]=size(X);
   for i=1:Nx     
     X(:,i)=(X(:,i)-para1(i))/para2(i);
   end
end
```

#### plsnipals

```
function [B,Wstar,T,P,Q,W,R2X,R2Y]=plsnipals(X,Y,A)
%+++ The NIPALS algorithm for both PLS-1 (a single y) and PLS-2 (multiple Y)
%+++ X: n x p matrix
%+++ Y: n x m matrix
%+++ A: number of latent variables
%+++ Code: Hongdong Li, lhdcsu@gmail.com, Feb, 2014
%+++ reference: Wold, S., M. Sj鰏tr鰉, and L. Eriksson, 2001. PLS-regression: a basic tool of chemometrics,
%               Chemometr. Intell. Lab. 58(2001)109-130.



varX=sum(sum(X.^2));
varY=sum(sum(Y.^2));
for i=1:A
    error=1;
    u=Y(:,1);
    niter=0;
    while (error>1e-8 && niter<1000)  % for convergence test
        w=X'*u/(u'*u);
        w=w/norm(w);
        t=X*w;
        q=Y'*t/(t'*t);  % regress Y against t;
        u1=Y*q/(q'*q);
        error=norm(u1-u)/norm(u);
        u=u1;
        niter=niter+1;
    end
    p=X'*t/(t'*t);
    X=X-t*p';
    Y=Y-t*q';
    
    %+++ store
    W(:,i)=w;
    T(:,i)=t;
    P(:,i)=p;
    Q(:,i)=q;
    
end

%+++ calculate explained variance
R2X=diag(T'*T*P'*P)/varX;
R2Y=diag(T'*T*Q'*Q)/varY;

Wstar=W*(P'*W)^(-1); 
B=Wstar*Q';
Q=Q';

%+++ 
```

#### plsda

```
function LDA=plslda(X,y,A,K,method,order,weight)
%+++ programmed according to NIPALS algorithm by Hongdong Li,Oct. 2006.
%+++ model: x=t*p'   y=t*r'=u*q'
%+++ y0 has to satisfy:  +1: positive class and -1:negative class;
%+++ A: number of latent variables
%+++ method:    'autoscaling', default
%        or 'center' 
%        or 'pareto';
%+++ weight: whether classes need to be weighted 
%            0: no weight, default
%            1: assigning the same weight to each class
%+++ Advisor: Yizeng Liang, yizeng_liang@263.net
%+++ H.D. Li, Feb. 8, 2009, lhdcsu@gmail.com
if nargin<7;weight=0;end
if nargin<6;order=1;end
if nargin<5;method='autoscaling';end;
if nargin<4;K=10;end
if nargin<3;A=2;end;
if nargin<2;warn('Wrong input parameters!');end
A=min([size(X) A]);
check=0; %+++ check data effectiveness;
if order==1
  [y,indexyy]=sort(y);
  X=X(indexyy,:);
elseif order==0
  indexyy=randperm(length(y));
  X=X(indexyy,:);
  y=y(indexyy);
end

[Mx,Nx]=size(X);%Mx行数146，Nx列数129
if K>Mx; K=Mx; end
A=min([size(X,1)-ceil(length(y)/K) size(X,2) A]);
yytest=[];YR=[];

groups = 1+rem(0:Mx-1,K);

yytest=[];YR=zeros(Mx,A);
B2=[];coef2=[];VIP2=[]; ws2=[];
for group=1:K
    testk = find(groups==group);  
    calk = find(groups~=group);
    Xcal=X(calk,:);ycal=y(calk);
    Xtest=X(testk,:);ytest=y(testk);
    %keyboard
    %data pretreatment    
    [Xcal,para1,para2]=pretreat(Xcal,method);
    %Xcal=Xcal(:,(sum(isnan(Xcal)))==0);
    Xcal(:,(sum(isnan(Xcal)))~=0)=0;
    ycals=ycal;  
    %ycals=pretreat(ycal,method);    
    Xtest=pretreat(Xtest,method,para1,para2);
    %Xtest=Xtest(:,(sum(isnan(Xtest)))==0);
    Xtest(:,(sum(isnan(Xtest)))~=0)=0;
%%%%%%%%%% PLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B,Wstar,T,P,Q,W,R2X,R2Y]=plsnipals(Xcal,ycals,A);
%[T2,c2,p2,u2,B2]=newpls(Xcal,ycals,A);

%%%%%%%%%% LDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:A
        %%%+++ train model.        
        TT=T(:,1:j); 
        C=ldapinv(TT,ycals,weight);                      %+++ Discriminant analysis
        %B0=[Wstar(:,1:j)*C(1:end-1);C(en d)];
        coef=[Wstar(:,1:j)*C(1:end-1);C(end)];  %=B0  
        coef2=[coef2 coef];
        VIP=vip(Xcal,ycals,T(:,1:j),W(:,1:j),Q(1:j,:));
        VIP2=[VIP2 VIP'];
        %+++ predict
       y_est=Xtest*coef(1:end-1)+coef(end);
        %Ttest=Xtest*Wstar(:,1:j);
        %y_est=Ttest*C(1:end-1)+C(end); 
       
        YR(testk,j)=y_est;
 end
  ws2=[ws2 Wstar];
end
  %+++ Original order
YR(indexyy,:)=YR;
y(indexyy)=y;
error=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:A
  error(i)=sum(sign(YR(:,i))~=y);    
end
error=error/Mx;

if sum(sum(isnan(T)))+sum(sum(isinf(T)))>0; check=1;end %+++ Check data effectivness;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if check==0;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mincv,index]=min(error);
ws=ws2(:,[index,index+A,index+2*A]);
VIP3=VIP2(:,[index,index+A,index+2*A]);
coef3=coef2(:,[index,index+A,index+2*A]);
coef=median(coef3,2);
VIP=median(VIP3,2);
wstar=median(ws,2);
B2=[B2 coef(1:end-1)];
ROC=roccurve(YR(:,index),y,0); %+++ ROC area.
%AUC=ROC.AUC;
k1=find(y==1);
k2=find(y==-1);
y_est=sign(YR(:,index));
sn=1-sum(y_est(k1)~=sign(y(k1)))/length(k1);
sp=1-sum(y_est(k2)~=sign(y(k2)))/length(k2);

%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
  LDA.method=method;
  %LDA.scale_para=[para1;para2];
  LDA.yreal=y;
  LDA.coef=coef;
  LDA.B=B2;
  %LDA.weight=W;
  LDA.wstar=wstar;
  %LDA.Xloadings=P;
  %LDA.R2X=R2X;
  %LDA.R2Y=R2Y;
  LDA.Sn=sn;
  LDA.Sp=sp;
  LDA.opt=index;
  %LDA.tpScores=tpt;
  %LDA.tpLoadings=tpp;
  %LDA.SR=SR;
  %LDA.check=0;
  %LDA.coef_lda_pc=C;
  %LDA.coef_lda_origin=[W(:,1:A)*C(1:end-1);C(end)];
  %LDA.yfit=yfit;
  LDA.VIP=VIP;
  LDA.y_est=y_est;
  LDA.error=1-ROC.accuracy;
  LDA.sensitivity=ROC.sensitivity;
  LDA.specificity=ROC.specificity;
  LDA.MCC=ROC.MCC;
  LDA.AUC=ROC.AUC;
elseif check==1
  LDA.check=1; 
end
%+++ There is a song you like to sing.
```

#### plotspread

```
function handles = plotSpread(varargin)
%PLOTSPREAD plots distributions of points by spreading them around the y-axis
%
% SYNOPSIS: handles = plotSpread(data,binWidth,spreadFcn,xNames,showMM,xValues)
%           handles = plotSpread(ah,...)
%
% INPUT data: cell array of distributions or nDatapoints-by-mDistributions
%           array, or array with data that is indexed by either
%           distributionIdx or categoryIdx, or both.
%       distributionIdx: grouping variable that determines to which
%           distribution a data point belongs. Grouping is
%           resolved by calling grp2idx, and unless xNames have
%           been supplied, group names determine the x-labels.
%           If the grouping variable is numeric, group labels also
%           determine x-values, unless the parameter xValues has
%           been specified.
%       distributionColors : color identifier (string, cell array of
%           strings), or colormap, with a single color, or one color per
%           distribution (or per entry in distributionIdx). Colors the
%           distributions. Default: 'b'
%       distributionMarkers : string, or cell array of strings, with either
%           a single marker or one marker per distribution (or per entry in
%           distributionIdx). See linespec for admissible markers.
%           Default: '.'
%		categoryIdx: grouping variable that determines group membership for data
%			points across distributions. Grouping is resolved by calling
%           grp2idx.
%       categoryColors : color identifier (cell array of
%           strings), or colormap, with one color per category.
%           Colors the categories, and will override distributionColors.
%           Default is generated using distinguishable_colors by Timothy E.
%           Holy.
%       categoryMarkers : cell array of strings, with one marker per
%           category. See linespec for admissible markers. Will override
%           distributionMarkers. Default: ''
%       binWidth : width of bins (along y) that control which data
%           points are considered close enough to be spread. Default: 0.1
%       spreadFcn : cell array of length 2 with {name,param}
%           if name is 'lin', the spread goes linear with the number of
%             points inside the bin, until it reaches the maximum of 0.9 at
%             n==param.
%           if name is 'xp', the spread increases as 1-exp(log(0.9)*x).
%             param is empty
%           Default {'xp',[]}
%       spreadWidth : width, along the x-axis (y-axis if flipped) that can
%           at most be covered by the points. Default:
%           median(diff(sort(xValues))); 1 if no xValues have been supplied
%       showMM : if 1, mean and median are shown as red crosses and
%                green squares, respectively. Default: 0
%                2: only mean
%                3: only median
%                4: mean +/- standard error of the mean (no median)
%                5: mean +/- standard deviation (no median)
%       xNames : cell array of length nDistributions containing x-tick names
%               (instead of the default '1,2,3')
%       xValues : list of x-values at which the data should
%                 be plotted. Default: 1,2,3...
%       xMode  : if 'auto', x-ticks are spaced automatically. If 'manual',
%                there is a tick for each distribution. If xNames is
%                provided as input, xMode is forced to 'manual'. Default:
%                'manual'.
%       xyOri  : orientation of axes. Either 'normal' (=default), or
%                'flipped'. If 'flipped', the x-and y-axes are switched, so
%                that violin plots are horizontal. Consequently,
%                axes-specific properties, such as 'yLabel' are applied to
%                the other axis.
%       yLabel : string with label for y-axis. Default : ''
%       ah  : handles of axes into which to plot
%
% OUTPUT handles: 3-by-1 cell array with handles to distributions,
%          mean/median etc, and the axes, respectively
%
% REMARKS: plotSpread is useful for distributions with a small number of
%          data points. For larger amounts of data, distributionPlot is
%          more suited.
%
% EXAMPLES: data = {randn(25,1),randn(100,1),randn(300,1)};
%           figure,plotSpread(data,[],[],{'25 pts','100 pts','300 pts'})
%
%            data = [randn(50,1);randn(50,1)+3.5]*[1 1];
%            catIdx = [ones(50,1);zeros(50,1);randi([0,1],[100,1])];
%            figure
%            plotSpread(data,'categoryIdx',catIdx,...
%                 'categoryMarkers',{'o','+'},'categoryColors',{'r','b'})
%
% END
%
% created with MATLAB ver.: 7.9.0.3470 (R2009b) on Mac OS X  Version: 10.5.7 Build: 9J61
%
% created by: jonas
% DATE: 11-Jul-2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def.binWidth = 0.1;
def.spreadFcn = {'xp',[]};
def.xNames = [];
def.showMM = false;
def.xValues = [];
def.distributionIdx = [];
def.distributionColors = 'b';
def.distributionMarkers = '.';
def.xMode = 'manual';
def.xyOri = 'normal';
def.categoryIdx = [];
def.categoryColors = [];
def.categoryMarkers = '';
def.yLabel = '';
def.spreadWidth = [];

%% CHECK INPUT

% check for axes handle
if ~iscell(varargin{1}) && length(varargin{1}) == 1 && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes')
    ah = varargin{1};
    data = varargin{2};
    varargin(1:2) = [];
    newAx = false;
else
    ah = gca;
    data = varargin{1};
    varargin(1) = [];
    % if the axes have children, it's not new (important for adjusting
    % limits below)
    newAx = isempty(get(ah,'Children'));
end

% optional arguments
parserObj = inputParser;
parserObj.FunctionName = 'plotSpread';
distributionIdx = [];distributionLabels = '';
if ~isempty(varargin) && ~ischar(varargin{1}) && ~isstruct(varargin{1})
    % old syntax
    parserObj.addOptional('binWidth',def.binWidth);
    parserObj.addOptional('spreadFcn',def.spreadFcn);
    parserObj.addOptional('xNames',def.xNames);
    parserObj.addOptional('showMM',def.showMM);
    parserObj.addOptional('xValues',def.xValues);
    
    parserObj.parse(varargin{:});
    opt = parserObj.Results;
    
    opt.distributionIdx = [];
    opt.distributionColors = def.distributionColors;
    opt.distributionMarkers = def.distributionMarkers;
    opt.xMode = def.xMode;
    opt.xyOri = def.xyOri;
    opt.categoryIdx = [];
    opt.categoryColors = def.distributionColors;
    opt.categoryMarkers = def.distributionMarkers;
    opt.yLabel = '';
    opt.spreadWidth = def.spreadWidth;
    
    for fn = fieldnames(def)'
        if ~isfield(opt,fn{1})
            % Manually adding the new defaults means a lot fewer bugs
            error('please add option %s to old syntax',fn{1});
        end
        if isempty(opt.(fn{1}))
            opt.(fn{1}) = def.(fn{1});
        end
    end
    
else
    % new syntax
    defNames = fieldnames(def);
    for dn = defNames(:)'
        parserObj.addParamValue(dn{1},def.(dn{1}));
    end
    
    
    parserObj.parse(varargin{:});
    opt = parserObj.Results;
end

% We want data to be a vector, so that indexing with both groupIdx and
% distributionIdx becomes straightforward, and so that we can conveniently
% eliminate NaNs that otherwise could mess up grouping.
% Consequently, if data is a cell array, we convert it, and build a
% corresponding distributionIdx (allowing a user-supplied distributionIdx
% to override, though), and then we go and take care of groupIdx. Once all
% three indices have been built, NaN can be removed.

if iscell(data)
    % make sure data is all n-by-1
    data = cellfun(@(x)x(:),data,'UniformOutput',false);
    nData = length(data);
    nn = cellfun(@numel,data);
    % make vector
    data = cat(1,data{:});
    distributionIdx = repeatEntries((1:nData)',nn);
else
    % distributions in columns
    nData = size(data,2);
    distributionIdx = repeatEntries((1:nData)',size(data,1));
    data = data(:);
end



% distribution groups
if ~isempty(opt.distributionIdx)
    [distributionIdx,distributionLabels,vals] = grp2idx(opt.distributionIdx);
    % convert data to cell array
    nData = length(distributionLabels);
    % if not otherwise provided, use group labels for xnames
    if isempty(opt.xNames)
        opt.xNames = distributionLabels;
        if ~iscell(opt.xNames)
            opt.xNames = num2cell(opt.xNames);
        end
    end
    if isnumeric(vals) && isempty(opt.xValues)
        opt.xValues = vals;
    end
end

if ~isempty(opt.xNames)
    opt.xMode = 'manual';
end


% distribution colors&markers
if ischar(opt.distributionColors)
    opt.distributionColors = {opt.distributionColors};
end
if iscell(opt.distributionColors)
    if length(opt.distributionColors) == 1
        % expand
        opt.distributionColors = repmat(opt.distributionColors,nData,1);
    elseif length(opt.distributionColors) ~= nData
        error('please submit one color per distribution (%i dist, %i colors)',nData,length(opt.distributionColors));
    end
    
else
    if size(opt.distributionColors,2) ~= 3
        error('please specify colormap with three columns')
    end
    if size(opt.distributionColors,1) == 1
        opt.distributionColors = repmat(opt.distributionColors,nData,1);
    elseif size(opt.distributionColors,1) ~= nData
        error('please submit one color per distribution (%i dist, %i colors)',nData,size(opt.distributionColors,1));
    end
    
    % create a cell array
    opt.distributionColors = mat2cell(opt.distributionColors,ones(nData,1),3);
end

if ischar(opt.distributionMarkers)
    opt.distributionMarkers = {opt.distributionMarkers};
end
if length(opt.distributionMarkers) == 1
    % expand
    opt.distributionMarkers = repmat(opt.distributionMarkers,nData,1);
elseif length(opt.distributionMarkers) ~= nData
    error('please submit one color per distribution (%i dist, %i colors)',nData,length(opt.distributionMarkers));
end


stdWidth = 1;
if isempty(opt.xValues)
    opt.xValues = 1:nData;
end

if isempty(opt.spreadWidth) 
    % scale width
    tmp = median(diff(sort(opt.xValues)));
    if ~isnan(tmp)
        stdWidth = tmp;
    end
else
    stdWidth = opt.spreadWidth;
end

if ~ischar(opt.xyOri) || ~any(ismember(opt.xyOri,{'normal','flipped'}))
    error('option xyOri must be either ''normal'' or ''flipped'' (is ''%s'')',opt.xyOri);
end


% check for categoryIdx/colors/markers
% If there are categories, check colors/markers individually first,
% then check whether any of them at all have been supplied, and
% if not, override distributionColors with default categoryColors

if isempty(opt.categoryIdx)
    categoryIdx = ones(size(distributionIdx));
    nCategories = 1;
    categoryLabels = '';
else
    [categoryIdx,categoryLabels] = grp2idx(opt.categoryIdx(:));
    nCategories = max(categoryIdx);
end

% plotColors, plotMarkers, plotLabels: nDist-by-nCat arrays
plotColors = repmat(opt.distributionColors(:),1,nCategories);
plotMarkers= repmat(opt.distributionMarkers(:),1,nCategories);

if isempty(distributionLabels)
    distributionLabels = opt.xNames;
    if isempty(distributionLabels)
        distributionLabels = cellstr(num2str(opt.xValues(:)));
    end
end

if nCategories == 1
    plotLabels = distributionLabels(:);
else
    plotLabels = cell(nData,nCategories);
    for iData = 1:nData
        for iCategory = 1:nCategories
            plotLabels{iData,iCategory} = ...
                sprintf('%s-%s',num2str(distributionLabels{iData}),...
                num2str(categoryLabels{iCategory}));
        end
    end
    
end




categoryIsLabeled = false;
if nCategories > 1
    % if not using defaults for categoryColors: apply them
    if ~any(strcmp('categoryColors',parserObj.UsingDefaults))
        if iscell(opt.categoryColors)
            if length(opt.categoryColors) ~= nCategories
                error('please supply one category color per category')
            end
            plotColors = repmat(opt.categoryColors(:)',nData,1);
            categoryIsLabeled = true;
        else
            if all(size(opt.categoryColors) ~= [nCategories,3])
                error('please supply a #-of-categories-by-3 color array')
            end
            plotColors = repmat( mat2cell(opt.categoryColors,ones(nCategories,1),3)', nData,1);
            categoryIsLabeled = true;
        end
    end
    
    if ~any(strcmp('categoryMarkers',parserObj.UsingDefaults))
        if length(opt.categoryMarkers) ~= nCategories
            error('please supply one category marker per category')
        end
        if ~iscell(opt.categoryMarkers)
            error('please supply a list of markers as cell array')
        end
        plotMarkers = repmat(opt.categoryMarkers(:)',nData,1);
        categoryIsLabeled = true;
    end
    
    if ~categoryIsLabeled
        % use distinguishable_colors to mark categories
        
        plotColors = repmat( mat2cell(...
            distinguishable_colors(nCategories),...
            ones(nCategories,1),3)', nData,1);
        
    end
    
end


% remove NaNs from data
badData = ~isfinite(data) | ~isfinite(distributionIdx) | ~isfinite(categoryIdx);
data(badData) = [];
distributionIdx(badData) = [];
categoryIdx(badData) = [];




%% TRANSFORM DATA
% Here, I try to estimate what the aspect ratio of the data is going to be
fh = figure('Visible','off');
if ~isempty(data)
    minMax = [min(data);max(data)];
else
    minMax = [0 1];
end
switch opt.xyOri
    case 'normal'
        plot([0.5;nData+0.5],minMax,'o');
    case 'flipped'
        plot(minMax,[0.5;nData+0.5],'o');
        
end
aspectRatio = get(gca,'DataAspectRatio');
close(fh);

tFact = aspectRatio(2)/aspectRatio(1);
if strcmp(opt.xyOri,'flipped')
    tFact = 1/tFact;
end

%% SPREAD POINTS
% assign either nData, or xValues number of values, in case we're working
% with group-indices
[m,md,sem,sd] = deal(nan(max(nData,length(opt.xValues)),1));
% augment data to make n-by-2
data(:,2) = 0;
for iData = 1:nData
    currentDataIdx = distributionIdx==iData;
    currentData = data(currentDataIdx,1);
    
    if ~isempty(currentData)
        
        % transform and sort
        currentData = currentData / tFact;
        %currentData = sort(currentData);
        
        % add x
        currentData = [ones(size(currentData))*opt.xValues(iData),currentData]; %#ok<AGROW>
        
        % step through the data in 0.1 increments. If there are multiple
        % entries, spread along x
        for y = min(currentData(:,2)):opt.binWidth:max(currentData(:,2))
            % find values
            valIdx = find(currentData(:,2) >= y & currentData(:,2) < y+opt.binWidth);
            nVal = length(valIdx);
            if nVal > 1
                % spread
                switch opt.spreadFcn{1}
                    case 'xp'
                        spreadWidth = stdWidth*0.9*(1-exp(log(0.9)*(nVal-1)));
                    case 'lin'
                        spreadWidth = stdWidth*0.9*min(nVal-1,opt.spreadFcn{2})/opt.spreadFcn{2};
                end
                spreadDist = spreadWidth / (nVal - 1);
                if isEven(nVal)
                    offset = spreadDist / 2;
                else
                    offset = eps;
                end
                for v = 1:nVal
                    currentData(valIdx(v),1) = opt.xValues(iData) + offset;
                    % update offset
                    offset = offset - sign(offset) * spreadDist * v;
                end
            end
        end
        
        % update data
        currentData(:,2) = data(currentDataIdx,1);
        data(currentDataIdx,:) = currentData;
        
        
        if opt.showMM > 0
            m(iData) = nanmean(currentData(:,2));
            md(iData) = nanmedian(currentData(:,2));
            sd(iData) = nanstd(currentData(:,2));
            sem(iData) = sd(iData)/sqrt(sum(isfinite(currentData(:,2))));
        end
    end % test isempty
end


%% plot
set(ah,'NextPlot','add')
ph = NaN(nData,nCategories);
for iData = 1:nData
    for iCategory = 1:nCategories
        currentIdx = distributionIdx == iData & categoryIdx == iCategory;
        if any(currentIdx)
            switch opt.xyOri
                case 'normal'
                    ph(iData,iCategory) = plot(ah,data(currentIdx,1),...
                        data(currentIdx,2),...
                        'marker',plotMarkers{iData,iCategory},...
                        'color',plotColors{iData,iCategory},...
                        'lineStyle','none',...
                        'DisplayName',plotLabels{iData,iCategory});
                case 'flipped'
                    ph(iData,iCategory) = plot(ah,data(currentIdx,2),...
                        data(currentIdx,1),...
                        'marker',plotMarkers{iData,iCategory},...
                        'color',plotColors{iData,iCategory},...
                        'lineStyle','none',...
                        'DisplayName',plotLabels{iData,iCategory});
            end
        end
    end
end


% if ~empty, use xNames
switch opt.xyOri
    case 'normal'
        switch opt.xMode
            case 'manual'
                set(ah,'XTick',opt.xValues);
                if ~isempty(opt.xNames)
                    set(ah,'XTickLabel',opt.xNames)
                end
            case 'auto'
                % no need to do anything
        end
        
        % have plot start/end properly
        minX = min(opt.xValues)-stdWidth;
        maxX = max(opt.xValues)+stdWidth;
        if ~newAx
            oldLim = xlim;
            minX = min(minX,oldLim(1));
            maxX = max(maxX,oldLim(2));
        end
        xlim([minX,maxX])
        
        ylabel(ah,opt.yLabel)
        
    case 'flipped'
        switch opt.xMode
            case 'manual'
                set(ah,'YTick',opt.xValues);
                if ~isempty(opt.xNames)
                    set(ah,'YTickLabel',opt.xNames)
                end
            case 'auto'
                % no need to do anything
        end
        
        % have plot start/end properly (for ease of copying, only switch
        % xlim to ylim
        minX = min(opt.xValues)-stdWidth;
        maxX = max(opt.xValues)+stdWidth;
        if ~newAx
            oldLim = ylim;
            minX = min(minX,oldLim(1));
            maxX = max(maxX,oldLim(2));
        end
        ylim([minX,maxX])
        
        xlabel(ah,opt.yLabel);
        
end


% add mean/median
mh = [];mdh=[];
if opt.showMM
    % plot mean, median. Mean is filled red circle, median is green square
    % I don't know of a very clever way to flip xy and keep everything
    % readable, thus it'll be copy-paste
    switch opt.xyOri
        case 'normal'
            if any(opt.showMM==[1,2])
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
            end
            if any(opt.showMM==[1,3])
                mdh = plot(ah,opt.xValues,md,'sg','MarkerSize',12);
            end
            if opt.showMM == 4
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,opt.xValues,m,sem);
            end
            if opt.showMM == 5
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,opt.xValues,m,sd);
            end
        case 'flipped'
            if any(opt.showMM==[1,2])
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
            end
            if any(opt.showMM==[1,3])
                mdh = plot(ah,md,opt.xValues,'sg','MarkerSize',12);
            end
            if opt.showMM == 4
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,m,opt.xValues,[sem,NaN(size(sem))]);
            end
            if opt.showMM == 5
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,m,opt.xValues,[sd,NaN(size(sd))]);
            end
    end
end

%==========================
%% CLEANUP & ASSIGN OUTPUT
%==========================

if nargout > 0
    handles{1} = ph;
    handles{2} = [mh;mdh];
    handles{3} = ah;
end
```

#### ldapinv

```
function C=ldapinv(X,y,weight)
%+++weight: =0   Bayesian approximation.
%+++        =1   Fisher DA, classes weighted.
%+++ The last element in C is the bias term.
%+++Hongdong Li,Nov. 23.

y=sign(y);
A=length(y);
B=length(find(y==1));
C=A-B;
r1=A/B;r2=A/C;R=[];
kp=find(y==1);
kn=find(y==-1);
XX=[[X(kp,:) ones(B,1)];-[X(kn,:) ones(C,1)]];

if weight==0
    BB=ones(A,1);
elseif weight==1
    BB=[ones(B,1)*r1;ones(C,1)*r2];
end

C=inv(XX'*XX)*XX'*BB;
```

#### isEven

```
function out = isEven(in)
%ISEVEN checks whether a number is even
%
% SYNOPSIS out = isEven(in)
%
% INPUT    in :  input (array) of numbers to be tested. 
% OUTPUT   out:  array of size(in) with 
%                   1 for even integers and zero
%                   0 for odd integers
%                 NaN for non-integers
%                out is a logical array as long as the input is all integers.
%
% c: jonas 5/05
% Last modified 11/24/2009 - Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = mod(in+1, 2);
% Set NaN for non-integer data, because they are neither odd or even
out((out ~= 0) & (out ~= 1)) = NaN;

% since doubles cannot be used for logical indexing, we should convert to
% logicals if possible. 
if all(isfinite(out(:)))
    out = logical(out);
end
```

