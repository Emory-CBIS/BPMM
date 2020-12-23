function [Graph,CP_all,Clusteri] = idPMAC(Yall,X,H,ncluster)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idPMAC.m: 
% BPMM(Bayesian Product Mxiture Modeling with covariates for precision matrix)
% Description: Fit a Bayesian product mixture model that imposes independent mixture priors 
% at each time scan and uses covariates to model the mixture weights, which results 
% in time-varying clusters of samples designed to pool information.
%
% Usage:
% [Graph,CP_all,sub_kmeans] = idPMAC(Yall,X,H,ncluster);
%
% Input:
% Yall          preprocessed fMRI data (T*V*Nsub), T is the total number of
%               scans, V is the number of ROIs, Nsub is the number of subjects
% X             covariates of all subjects (Nsub*C), C is the number of
%               covariates for subjects 
% H             number of mixture compoment (default value is round(nsub/8))
% ncluster      number of clusters, default value is same as H
%
% Output: 
% Graph         the estiamted pairwise correlation (V*V*T*Nsub)
% CP_all        change points for all subjects
% Clusteri      estimated cluster of all subjects
% 
% AUTHORS:
% Suprateek Kundu     Emory University
% Jin Ming            Emory University
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data load
if nargin < 2
    error('more input arguments needed.');
end

nsub=size(Yall,3);
if nargin < 3 || isempty(H)
    H = round(nsub/8);
end

if nargin < 4 || isempty(ncluster)
   ncluster = H;
end

fprintf('running idPMAC \n');
dx = size(X,2);
T=size(Yall,1);
V=size(Yall,2);
nsub=size(Yall,3);
dis_min = 20;
kappa = zeros(V,V,T,nsub);
for k=1:nsub
    Y = Yall(:,:,k);
    for i=(dis_min+1):(T-dis_min)
    Yi = Y((i-dis_min):(i+dis_min),:);
    kappa(:,:,i,k) = corr(Yi);
    end
    for i=1:dis_min
kappa(:,:,i,k) = kappa(:,:,dis_min+1,k);
    end
    for i=(T-dis_min+1):T
kappa(:,:,i,k) = kappa(:,:,T-dis_min,k);
    end
kappa(1:(V+1):end) = 0;
end

vary2 = zeros(nsub,V);
for k=1:nsub
for i=1:V
    Y=Yall(:,:,k);
   vary2(k,i) = var(Y(:,i)); 
end
end
vary = mean(vary2,2)';
varysub = zeros(V,V,T,nsub);
for i=1:(V-1)
    for j=(i+1):V
        for t=1:T
  varysub(i,j,t,:)=vary;
        end
    end
end
Smatrix = zeros(V,V,T,nsub);
Omega = zeros(V,V,T,nsub);
r1 = 0.01;
for sub=1:nsub
   for t=1:T
       Yi = squeeze(Yall(t,:,sub));
   Smatrix(:,:,t,sub)=Yi'*Yi;
   s = Yi'*Yi;
  [Xi,W,~,~,~,~] = QUIC('default', s, r1, 1e-6, 2, 100);
  Omega(:,:,t,sub) = Xi;
   end
end

OmegaH = zeros(V,V,T,H);

xiH = zeros(V,T,nsub,H);
psiH = zeros(V,T,nsub,H);
betaH = zeros(V,T,dx+1,H);
for i=1:V
   for t=1:T
      omegahi = squeeze(Omega(i,:,t,:));
      [ind1,ind] = sort(omegahi(1,:));
      ind2 = quantile(1:nsub,H-1);
      ind2 = round([0,ind2,nsub]);
      for j=1:H
         indi1 = ind2(j)+1;
         indi2 = ind2(j+1);
         indi3 = ind(indi1:indi2);
         OmegaH(i,:,t,j)=mean(omegahi(:,indi3),2);
      end
   end
end

sigmaH = ones(1,H);
for h=1:H
   sigmaH(h) = std(reshape(OmegaH(:,:,:,h),1,[])) ;
end

for i=1:V
    Bsave = ones(dx+1,4);
    for t=1:T
        xi = ones(nsub,H)*(1/H); %subject-based
        xiH(i,t,:,:) = xi;
    betaH(i,t,:,1:(H-1)) = Bsave;
    B = squeeze(betaH(i,t,:,:));
    psi = exp(X1*B)./(sum(exp(X1*B),2)); 
    while sum(sum(isnan(psi)))>0
       Bnew = B/10; 
       B = Bnew;
       psi = exp(X1*B)./(sum(exp(X1*B),2)); %subject-average
    end
     psiH(i,t,:,:) = psi;     
    end
end
alpha = 1;
Kappa = zeros(V,T,nsub);
for i=1:V
    for t=1:T
        for sub=1:nsub
   omegai = Omega(i,i,t,sub);
   omegaj = Omega(i,:,t,sub);
   omegaj(i) = [];
   OmegaJ = Omega(:,:,t,sub);
   OmegaJ(i,:) = [];
   OmegaJ(:,i) = [];
   Kappa(i,t,sub) = omegai - omegaj*inv(OmegaJ)*omegaj';   
        end
    end
end

lambda = 1;
sub_con = zeros(nsub,nsub,T);
Hsave = zeros(V,V,T);
dis_min = 20;
maxits = 40;
niters = 1;
 a1 = 1;
 b1 = 1/10;
sigma = 100;
sigmaB = diag(ones(dx+1,1)*1000);

gammabar = zeros(V,V,T,H);
whjlt = zeros(V,V,T,H);
%target function
obj_check = zeros(maxits+1,8);

objfun = zeros(maxits+1,1);
err=zeros(maxits+1,1);
err(1,1)=1;
eps = 6e-4;
obj1 = 0;
for sub=1:nsub
    [sub]
   for t=1:T
       xiHi = xiH(:,t,sub,:);
      yit = squeeze(Yall(t,:,sub)); 
      omeyit = Omega(:,:,t,sub);
      if det(omeyit) > 0
      pyit1 = log(det(omeyit));
      else 
          pyit1 = 0;
      end
      pyit2 = 0;
      pyit3 = 0;
      pyit4 = 0;
      smatrixi = Smatrix(:,:,t,sub);
      for i=1:V
      xiHi = squeeze(xiH(i,t,sub,:));
      omevit = omeyit(i,:);
      omevit(i) = [];
      Omevit = omeyit;
      Omevit(i,:) = [];
      Omevit(:,i) = [];
      kvit = omeyit(i,i) - omevit * Omevit * omevit';
      if kvit>0
      pyit2 = pyit2 + log(kvit);
      end
      svit = smatrixi(i,:);
      svit(i) = [];
      pyit3a = ((smatrixi(i,i)+alpha)/2)*kvit;
      pyit3b = 0.5* omevit*(1*eye(V-1) + (smatrixi(i,i)+alpha)*Omevit)*(omevit');
      pyit3c = 2*svit*omevit';
      pyit3 = pyit3 - pyit3a - pyit3b + pyit3c;
        for h=1:H
            sigmah = sigmaH(h);
      mustarhi = OmegaH(i,:,t,h);
      mustarhi(i) = [];
      omevitj = omeyit(i,:);
      omevitj(i) = [];
      pyit4 = pyit4 + (xiHi(h)/sigmah) * (omevitj-mustarhi)*(omevitj-mustarhi)';
        end
      end
      obj1 = obj1+pyit1+pyit2+pyit3+pyit4;
   end
end

obj2 = 0;
for h=1:H
    sigmah = sigma;
    xi2 = xiH(:,:,:,h);
    psi2 = psiH(:,:,:,h);
    obj2 = obj2 + sum(sum(nansum(xi2)))*log(sigmah)+ sum(sum(nansum(xi2.* log(psi2).*(psi2>0))));
end


obj3 = lambda*sum(sum(sum(nansum(OmegaH(:,:,2:T,:)-OmegaH(:,:,1:T-1,:)))));
obj4 = 0;
for h=1:H
    sigmah = sigmaH(h);
obj4 = obj4 + (a1+1)*log(sigmah) + b1/sigmah;
end
objfun(niters) = obj1 + obj2 + obj3 + obj4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter updates
while niters < maxits
    niters = niters+1;
    sub_con = zeros(nsub,nsub,T);   
    fprintf('Parameters update %2.0f \n',niters)

for i=1:V
        for t=1:T
     H =5;
    omegahi = squeeze(Omega(i,:,t,:));
    omegahj =squeeze(OmegaH(i,:,t,:));
    xi = squeeze(xiH(i,t,:,:)); %subject-based
    psi = squeeze(psiH(i,t,:,:));   %subject-average
    xiit = 0;
         %E-step
         xiit = xiit +1;   
         xinormpdf = zeros(nsub,H);
        for h=1:H
            sigmah = sigmaH(h)*T;
            xinormpdf(:,h) = mvnpdf(omegahi',omegahj(:,h)',sqrt(sigmah)*eye(V));
        end
        nxi = xi;
        for k=1:nsub
        nxi(k,:) = psi(k,:).*xinormpdf(k,:)/(sum(psi(k,:).*xinormpdf(k,:)));
        end
    if sum(sum(isnan(nxi)))==0
       xi=nxi; 
    end
    xi = xi./(sum(xi,2));
    Hsave(i,j,t) = H;
    xiH(i,t,:,:)=xi;
     Q = xi*xi';
     Q(isnan(Q)) = 0;
    sub_con(:,:,t) = sub_con(:,:,t) +Q;   
        end
        for h=1:(H-1)
           zih = zeros(nsub,T);
           wih = zeros(nsub,T);
           for t=1:T
               nxi = squeeze(xiH(i,t,:,:));
         B = squeeze(betaH(i,t,:,:));
        npsi = exp(X1*B)./(sum(exp(X1*B),2));
         while sum(sum(isnan(npsi)))>0
       Bnew = B/10; 
       B = Bnew;
       npsi = exp(X1*B)./(sum(exp(X1*B),2)); 
          end
        probh = npsi(:,h);
        probh(probh==1) = 1-1e-5;
        probh(probh==0) = 1e-5;
        bbb = X1*B;
        zih(:,t) = bbb(:,h) + (nxi(:,h)-probh)./(probh.*(1-probh));
        wih(:,t) = probh.*(1-probh);
           end
              A1 = zeros(dx+1,dx+1);
        A2 = zeros(dx+1,1);
            for sub=1:nsub
            wihi = sum(wih(sub,:));
            xixit = X1(sub,:)' * X1(sub,:);
            A1 = A1+ wihi* xixit;
            wihzih = nansum(wih(sub,:).*zih(sub,:));
            A2 = A2 + wihzih * X1(sub,:)';
            end
        A1 = A1 + inv(sigmaB);
        betai  = inv(A1)*A2;
            for t=1:T
       betaH(i,t,:,h) = betai;
            end
        end
        B = squeeze(betaH(i,t,:,:));
    npsi = exp(X1*B)./(sum(exp(X1*B),2));
    while sum(sum(isnan(npsi)))>0
       Bnew = B/10; 
       B = Bnew;
       npsi = exp(X1*B)./(sum(exp(X1*B),2)); 
    end
    for t=1:T
    psiH(i,t,:,:) = npsi;
    end
end   %end updates of Omega

%update Omega
for sub=1:nsub
    for t=1:T
        burnin  = 1; nmc = 2;
a_lambda = 1; b_lambda = 0.1; % Prior hyperparameters
S = Smatrix(:,:,t,sub);
C = Omega(:,:,t,sub); % Initial values 
Sig = inv(C);
sumD =sum(squeeze(xiH(:,t,sub,:)),2);
sumD2 = zeros(V-1,V);
for ii=1:V
    xiHi = squeeze(xiH(:,t,sub,:));
    xiHii = xiHi(ii,:);
    omegaHii=squeeze(OmegaH(ii,:,t,:));
    omegaHii(ii,:)=[];
    sumD2(:,ii)=sum(xiHii.*omegaHii,2);
end
[Sig_save,C_save,lambda_save] = BayesGLasso_Columnwise_Jin(S,1,Sig,C,sumD,sumD2,a_lambda,b_lambda,burnin,nmc);
Omega(:,:,t,sub)=squeeze(mean(Sig_save,3));
    end
end

for h=1:5
    sigmah = sigmaH(h);
    for i=1:(V-1)
        for j=(i+1):V
            for t=1:T
whjlt(i,j,t,h) =  sqrt(sum(xiH(i,t,:,h))/(2*sigmah));
whjlt(j,i,t,h) = whjlt(i,j,t,h);
gammabar(i,j,t,h) = whjlt(i,j,t,h)*sum(squeeze(Omega(i,j,t,:)).*(squeeze(xiH(i,t,:,h))))/(sum(xiH(i,t,:,h)));
gammabar(j,i,t,h) = gammabar(i,j,t,h);
            end
        end
    end
end
gammabar(isnan(gammabar)) = 0;
Mx = zeros(T,T);

%update of gammaH
for h=1:H
    mu2 = zeros(V,V,T);
for i=1:V
    ind3 = 1:V;
    ind3(i) = [];
     for j=ind3
    weight2 = squeeze((whjlt(i,j,:,h)));
    gamma2 = squeeze(gammabar(i,j,:,h));
    for t=1:T
        Mx(t,1:t) = weight2(t);
    end
    for tt=1:T
       if isnan(gamma2(tt))== 1
        gamma2(tt)=(gamma2(tt-1)+gamma2(tt+1))/2;
        if isnan(gamma2(tt))==1
            gamma2(tt)=gamma2(tt-1);
        end
       end
    end
            [B] = glmnet(Mx,gamma2);
            nlambda = size(B.beta,2);
            BICforlasso = zeros(1,nlambda);
            for kk=1:nlambda
                B1 = B.beta(:,kk);
                B0 = B.a0(kk);
                bdf = B.df(kk) +1 ;
                lambdaforlasso1 = B.lambda(kk);
                BICforlasso(kk) = bdf*log(T) + 2*(sum(((gamma2 - Mx*B1-B0).^2))) + 2*lambdaforlasso1*sum(abs(B1));
            end
              [~,idxLambda1SE] = min(BICforlasso);
                if length(idxLambda1SE)~=1
                idxLambda1SE = min(idxLambda1SE); 
                end
                B2 = B.beta(:,idxLambda1SE);
                B0 = B.a0(idxLambda1SE);
                B2(1) = B2(1) + B0/mean(Mx(:,1));
                   res = cumsum(B2);
OmegaH(i,j,:,h)=res; 
OmegaH(j,i,:,h)=res;
        end
    end
end    

for h=1:H
   sigmaH(h) = std(reshape(OmegaH(:,:,:,h),1,[])) ;
end

siggamma = 0.1;
obj1 = 0;
for sub=1:nsub
   for t=1:T
       xiHi = xiH(:,t,sub,:);
      yit = squeeze(Yall(t,:,sub)); 
      omeyit = Omega(:,:,t,sub);
      if det(omeyit) > 0
      pyit1 = log(det(omeyit));
      else 
          pyit1 = 0;
      end
      pyit2 = 0;
      pyit3 = 0;
      pyit4 = 0;
      smatrixi = Smatrix(:,:,t,sub);
      for i=1:V
      xiHi = squeeze(xiH(i,t,sub,:));
      omevit = omeyit(i,:);
      omevit(i) = [];
      Omevit = omeyit;
      Omevit(i,:) = [];
      Omevit(:,i) = [];
      kvit = omeyit(i,i) - omevit * Omevit * omevit';
      if kvit>0
      pyit2 = pyit2 + log(kvit);
      end
      svit = smatrixi(i,:);
      svit(i) = [];
      pyit3a = ((smatrixi(i,i)+alpha)/2)*kvit;
      pyit3b = 0.5* omevit*(1*eye(V-1) + (smatrixi(i,i)+alpha)*Omevit)*(omevit');
      pyit3c = 2*svit*omevit';
      pyit3 = pyit3 - pyit3a - pyit3b + pyit3c;
        for h=1:H
            sigmah = sigmaH(h);
      mustarhi = OmegaH(i,:,t,h);
      mustarhi(i) = [];
      omevitj = omeyit(i,:);
      omevitj(i) = [];
      pyit4 = pyit4 + (xiHi(h)/sigmah) * (omevitj-mustarhi)*(omevitj-mustarhi)';
        end
      end
      obj1 = obj1+pyit1+pyit2+pyit3+pyit4;
   end
end
obj2 = 0;
for h=1:H
    sigmah = sigma;
    xi2 = xiH(:,:,:,h);
    psi2 = psiH(:,:,:,h);
    obj2 = obj2 + sum(sum(nansum(xi2)))*log(sigmah)+ sum(sum(nansum(xi2.* log(psi2).*(psi2>0))));
end
obj3 = lambda*sum(sum(sum(nansum(OmegaH(:,:,2:T,:)-OmegaH(:,:,1:T-1,:)))));
obj4 = 0;
for h=1:H
    sigmah = sigmaH(h);
obj4 = obj4 + (a1+1)*log(sigmah) + b1/sigmah;
end
objfun(niters) = obj1 + obj2 + obj3 + obj4;

% if isreal(objfun(niters))==0
%    fprintf('Error Non-real objfun \n') 
% end

err(niters)=(objfun(niters-1)-objfun(niters))/abs(objfun(niters));
if abs(err(niters)) < eps && err(niters)<0 && niters >25
    break
end
end


%% subject clustering
sub_kmeans = zeros(nsub,T);
%save('withcoef1.mat','sub_con','gamma','SigTrue','gammaH','xiH','psiH','betaH')
%ss
fprintf('Subject culstering \n');
for t=1:T
sub_coni = sub_con(:,:,t);
sub_grou = kmeans(sub_coni,ncluster);
sub_grou2 = sub_grou;
sub_id = 1:nsub;
groum = zeros(1,ncluster);
for tt=1:ncluster
   groum(tt)=mean(sub_id(sub_grou==tt)) ;
end
[~,gind] = sort(groum);
for tt=1:ncluster
sub_grou2(sub_grou==gind(tt))=tt;
end
sub_kmeans(:,t)=sub_grou2;
end

samegroupp = zeros(nsub,nsub);
for i=1:(nsub-1)
    for j=(i+1):nsub
    samegroupp(i,j)= sum(sub_kmeans(i,:)==sub_kmeans(j,:))/T;
    samegroupp(j,i)=samegroupp(i,j);
    end
end
Clusteri = kmeans(samegroupp,ncluster);

Graph = Omega;
%% Network CP
dy_cp = zeros(nsub,T);
dy_cp2 = zeros(V,V,T,nsub);
for sub=1:nsub
    [sub]
   for i=1:(V-1)
       for j=(i+1):V
        ploti = squeeze(Omega(i,j,:,sub))';    %change the input T data
        res2 = simpleGFL2(ploti',weight3); 
        est_edgei = res2.jumps;
        for t=1:T
              if ismember(t,est_edgei)==1
                    dy_cp(sub,t) = dy_cp(sub,t) + 1;
                    dy_cp2(i,j,t,sub) = 1;
              end
        end
       end
    end
end

perceni = 0.5;
maxdist = 2;
CP_all = zeros(nsub,20);
for sub=1:nsub
    cp_esti = dy_cp(sub,:);
    cp_esti(T)=0;
    cp_est = cp_clus(cp_esti,perceni,maxdisti);
    cp_est(cp_est==1) = [];
    cp_est(cp_est==T) = [];
    CP_all(sub,1:length(cp_est)) = cp_est;
end

end
