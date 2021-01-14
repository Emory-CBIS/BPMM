function [Graph,CP_all,Clusteri] = BPMMPR(Yall,H,ncluster)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BPMMC_PR.m: 
% BPMM-PC (Bayesian Product Mxiture Modeling without covariates for precision matrix)
% Description: Fit a Bayesian product mixture model that imposes independent mixture priors 
% at each time scan, which results in time-varying clusters of samples designed to pool information.
%
% Usage:
% [Graph,CP_all,sub_kmeans] = BPMM_PR(Yall,H);
%
% Input:
% Yall          preprocessed fMRI data (T*V*Nsub), T is the total number of
%               scans, V is the number of ROIs, Nsub is the number of subjects
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
if nargin < 1
    error('more input arguments needed.');
end

nsub=size(Yall,3);
if nargin < 2 || isempty(H)
    H = round(nsub/8);
end

if nargin < 3 || isempty(ncluster)
   ncluster = H;
end

fprintf('running BPMM-PR \n');

T=size(Yall,1);
V=size(Yall,2);
dis_min = 20;

a=[1:T-1]';
weight3 = sqrt(T./(a.*(T-a)));
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
for t=1:T
    for i=1:(V-1)
        for j=(i+1):V
        omegahi = squeeze(Omega(i,j,t,:));
    omegahj =quantile(omegahi,H);
    if length(unique(omegahj)) == 1
        aa = [-0.2, -0.1, 0 , 0.1 , 0.2];
        omegahj = omegahj + aa;
    end
    OmegaH(i,j,t,:) = omegahj;
        xi = ones(nsub,H)*(1/H); %subject-based
        xiH(i,t,:,:) = xi;
    psi = ones(1,H)*(1/H);   %subject-average
     for sub=1:nsub
     psiH(i,t,sub,:) = psi;
     end
        end
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
   Kappa(i,t,sub) = omegai - omegaj*OmegaJ*omegaj';   
        end
    end
end

lambda = 1;
sub_con = zeros(nsub,nsub,T);
Hsave = zeros(V,V,T);

maxits = 100;
niters = 1;
 a1 = 1;
 b1 = 1/10;
sigma = 1;

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
      pyit3b = 0.5* omevit*(sigma*eye(V-1) + (smatrixi(i,i)+alpha)*Omevit)*(omevit');
      pyit3c = 2*svit*omevit';
      pyit3 = pyit3 - pyit3a - pyit3b + pyit3c;
        for h=1:H
      mustarhi = OmegaH(i,:,t,h);
      mustarhi(i) = [];
      omevitj = omeyit(i,:);
      omevitj(i) = [];
      pyit4 = pyit4 + (xiHi(h)/sigma) * (omevitj-mustarhi)*(omevitj-mustarhi)';
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
obj4 = (a1+1)*log(sigma) + b1/sigma;
objfun(niters) = obj1 + obj2 + obj3 + obj4;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameter updates


while niters < maxits
    niters = niters+1;
    sub_con = zeros(nsub,nsub,T);   
    fprintf('Parameters update %2.0f \n',niters)

for i=1:V
        for t=1:T
    omegahi = squeeze(Omega(i,:,t,:));
    omegahj =squeeze(OmegaH(i,:,t,:));
    xi = squeeze(xiH(i,t,:,:)); %subject-based
    psi = squeeze(psiH(i,t,1,:))';   %subject-average
    xiit = 0;
    sigmahi = sigma;
        %E-step
         xiit = xiit +1;   
         xinormpdf = zeros(nsub,H);
        for h=1:H
            xinormpdf(:,h) = mvnpdf(omegahi',omegahj(:,h)',sqrt(sigmahi)*eye(V));
        end
        nxi = xi;
        for k=1:nsub
        nxi(k,:) = psi.*xinormpdf(k,:)/(sum(psi.*xinormpdf(k,:)));
        end
     npsi = nanmean(nxi,1);
     if (sum(isnan(npsi)>0))
        npsi = ones(1,H)*(1/H);
     end
     if sum(npsi)==0
        npsi =  ones(1,H)*(1/H);
     end
     if sum(sum(isnan(nxi)))==0
        xi = nxi;
     end
    Hsave(i,j,t) = H;
    xiH(i,t,:,:)=xi;
    psiH(i,t,:,:)=repmat(npsi,nsub,1);
    
    %update of Omega
    for sub=1:nsub
        sigmaomegah = eye(V-1);
        svvi = Smatrix(i,i,t,sub);
        omegai1 = Omega(:,:,t,sub);
        omegai1(i,:) = [];
        omegai1(:,i) = [];
        svt = Smatrix(i,:,t,sub);
        svt(i) = [];
        omeright = 2*svt;
        right2 = zeros(1,V-1);
        left2 = 0;
        for h=1:H
            xij = xi(sub,h);
            omegahj = squeeze(OmegaH(i,:,t,h));
            omegahj(i) = [];
           right2 =  xij * omegahj + right2;
           left2 = left2 + xij;
        end
        left2 = left2*eye(V-1);
        omeupdat = inv(sigmaomegah+(svvi+alpha)*omegai1+left2) * ((2*svt+right2)');
        ind = 1:V;
        ind(i) = [];
        Omega(i,ind,t,sub) = omeupdat;
        Omega(ind,i,t,sub) = omeupdat;
    end
     
     Q = xi*xi';
     Q(isnan(Q)) = 0;
    sub_con(:,:,t) = sub_con(:,:,t) +Q;   
        end
        %update Kappa
for sub=1:nsub
   svvt = Smatrix(i,i,t,sub);
   nkappai = 2/(svvt+alpha);
   Kappa(i,t,sub) = nkappai;
   omegai1 = Omega(:,:,t,sub);
        omegai1(i,:) = [];
        omegai1(:,i) = [];
        omega2l = Omega(i,:,t,sub);
        omega2l(i) = [];
   omediagi = nkappai + omega2l*omegai1*omega2l';
   Omega(i,i,t,sub) = omediagi;
end

end   %end updates of Omega


for h=1:H
    for i=1:(V-1)
        for j=(i+1):V
            for t=1:T
whjlt(i,j,t,h) =  sqrt(sum(xiH(i,t,:,h))/(2*sigma));
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
      pyit3b = 0.5* omevit*(sigma*eye(V-1) + (smatrixi(i,i)+alpha)*Omevit)*(omevit');
      pyit3c = 2*svit*omevit';
      pyit3 = pyit3 - pyit3a - pyit3b + pyit3c;
        for h=1:H
      mustarhi = OmegaH(i,:,t,h);
      mustarhi(i) = [];
      omevitj = omeyit(i,:);
      omevitj(i) = [];
      pyit4 = pyit4 - 0.5* (xiHi(h)/sigma) * (omevitj-mustarhi)*(omevitj-mustarhi)';
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
    obj2 = obj2 - sum(sum(nansum(xi2)))*log(sigmah)+ sum(sum(nansum(xi2.* log(psi2).*(psi2>0))));
end

obj3 = lambda*sum(sum(sum(nansum(OmegaH(:,:,2:T,:)-OmegaH(:,:,1:T-1,:)))));
obj4 = (a1+1)*log(sigma) + b1/sigma;
objfun(niters) = obj1 + obj2 + obj3 + obj4;
      
if isreal(objfun(niters))==0
   fprintf('Error Non-real objfun \n') 
end

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
