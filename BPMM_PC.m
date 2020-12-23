function [Graph,CP_all,Clusteri] = BPMMPC(Yall,H,ncluster)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BPMMC_PC.m: 
% BPMM-PC (Bayesian Product Mxiture Modeling without covariates for pairwise correlation)
% Description: Fit a Bayesian product mixture model that imposes independent mixture priors 
% at each time scan, which results in time-varying clusters of samples designed to pool information.
%
% Usage:
% [Graph,CP_all,sub_kmeans] = BPMMPC(Yall,H);
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

fprintf('running BPMM-PC \n');

%% parameters set-up
T=size(Yall,1);
V=size(Yall,2);
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
gamma = 0.5 * log((1+kappa)./(1-kappa));

for i=2:V
    for j=1:(i-1)
            gamma(i,j,:,:) = 0;
    end
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
        
Hsave = zeros(V,V,T);

maxits = 100;
niters = 1;
a1 = 1;
b1 = 1/10;

eta = 0;
dimd = 4;  
Sigma_u = diag(ones(dimd,1));
mu_u = zeros(dimd,1);
uk = zeros(dimd,V);
for i=1:V
    uk(:,i)=mvnrnd(mu_u,Sigma_u)';
end
ukH = zeros(dimd,V,H);
for h=1:H
   ukH(:,:,h) = uk; 
end
Lambda = diag(ones(dimd,1));
tau2 = 10;
nu = zeros(V,V,T);
piT = zeros(V,V,T);
for i=1:(V-1)
    for j=(i+1):V
        for t=1:T
    meani  = uk(:,i)'*Lambda*(uk(:,j));
    nui = mvnrnd(meani,tau2);
    nu(i,j,t)= nui;  
    pii = exp(nui)/(1+exp(nui));
    piT(i,j,t)=pii;
        end
    end
end
nuH = zeros(V,V,T,H);
piTH = zeros(V,V,T,H);
for h=1:H
   nuH(:,:,:,h) = nu;
   piTH(:,:,:,h)= piT;
end
zH = binornd(1,piTH);

% initial guess of the Gamma (Fisher-transfered pairwise corr)
gammaH = zeros(V,V,T,H);
xiH = zeros(V,V,T,nsub,H);
psiH = zeros(V,V,T,nsub,H);
for i=1:(V-1)
    for j=(i+1):V
        for t=1:T
    gammahi = squeeze(gamma(i,j,t,:));
    gammahj =quantile(gammahi,H);
    xi = ones(nsub,H)*(1/H); %subject-based
    psi = ones(1,H)*(1/H);   %subject-average
    gammaH(i,j,t,:) = gammahj;
    xiH(i,j,t,:,:) = xi;
    for sub=1:nsub
    psiH(i,j,t,sub,:) = psi;
    end
        end
    end
end
gammabar = zeros(V,V,T,H);
whjlt = zeros(V,V,T,H);
siggamma = 0.1;
for h=1:5
    for i=1:(V-1)
        for j=(i+1):V
            for t=1:T
whjlt(i,j,t,h) =  sqrt(sum(xiH(i,j,t,:,h))/(2*siggamma));
gammabar(i,j,t,h) = whjlt(i,j,t,h)*sum(gamma(i,j,t,:).*xiH(i,j,t,:,h))/(sum(xiH(i,j,t,:,h)));
            end
        end
    end
end
% target function
objfun = zeros(maxits+1,1);
err=zeros(maxits+1,1);
err(1,1)=1;
eps = 1e-4;
yit2 = zeros(V,V,T,nsub);
yitjt = zeros(V,V,T,nsub);
for sub=1:nsub
    for t=1:T
   Yi = Yall(t,:,sub);
   yitjt(:,:,t,sub) = Yi'*Yi;
   yit2(:,:,t,sub)=Yi.^2 + Yi.^2';
    end
end
obj1 = -0.5*log(1-((exp(2*gamma)-1)./(exp(2*gamma)+1)).^2)-(yit2 - 2*((exp(2*gamma)-1)./(exp(2*gamma)+1)).*yitjt)./(2*varysub.*(1-((exp(2*gamma)-1)./(exp(2*gamma)+1)).^2));
obj1(obj1==-Inf) = 0;
obj2 = sum(sum(sum(nansum(obj1))));
obj4 = 0;
for h=1:H
   xi2 = xiH(:,:,:,:,h);
   psi2 = psiH(:,:,:,:,h);
   obj3 = 0.5*(xi2.*(gamma-gammaH(:,:,:,h)).^2)/(siggamma)+ 0.5*gamma*log(siggamma)+ xi2.* log(psi2); %- xi2.*log(1-psi2);
   obj4 = obj4 + sum(sum(sum(nansum(obj3))));
end
obj5 = 0;
lambda = 1;
obj5 =obj5 + lambda*sum(sum(sum(nansum(gammaH(:,:,2:T,:)-gammaH(:,:,1:(T-1),:))))) - (a1+1)*log(siggamma)-b1/siggamma ;

objfun(niters) = obj2 + obj4 + obj5;
subgroup = zeros(V,V,T,nsub);           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate of Gamma and GammaH
% parameter updates

while niters < maxits
    niters = niters+1;
    sub_con = zeros(nsub,nsub,T);   
    fprintf('Parameters update %2.0f \n',niters)
for i=1:(V-1)
    for j=(i+1):V
        for t=1:T
     H =5;
    gammahi = squeeze(gamma(i,j,t,:));
    gammahj =squeeze(gammaH(i,j,t,:));
    xi = squeeze(xiH(i,j,t,:,:)); %subject-based
    psi = squeeze(psiH(i,j,t,1,:))';   %subject-average
    xiit = 0;
    sigmahi = siggamma;
        %E-step
        %Update xi and psi
         xinormpdf = zeros(nsub,H);
        for h=1:H
            xinormpdf(:,h) = normpdf(gammahi,gammahj(h),sqrt(sigmahi));
        end
        nxi = xi;
        for k=1:nsub
        nxi(k,:) = psi.*xinormpdf(k,:)/(sum(psi.*xinormpdf(k,:)));
        end
         if sum(sum(isnan(nxi))) > 0
            [i j t]
             break
         end
     npsi = nanmean(nxi,1);
     if (sum(isnan(npsi)>0))
        npsi = ones(1,H)*(1/H);
     end
        xi = nxi;
    Hsave(i,j,t) = H;
    xiH(i,j,t,:,:)=xi;
    psiH(i,j,t,:,:)=repmat(npsi,nsub,1);
       xiit2 = 1;xidiff2 =1; lim1 =1e-2;
       %Update gamma 
       while xiit2 < 100 && xidiff2 > lim1
     xiit2 = xiit2 +1;
     yt2 = zeros(nsub,1);
     ytp = zeros(nsub,1); 
     h1 = zeros(nsub,1);
     for kk1=1:nsub
         yt2(kk1) = Yall(t,i,kk1)^2 + Yall(t,j,kk1)^2;
         ytp(kk1) = Yall(t,i,kk1)*Yall(t,j,kk1);
         h1(kk1) = sum(xi(kk1,:)'.*(gammahi(kk1)-gammahj)./sigmahi);
     end
     h2 = sum(xi./sigmahi,2);
     d1 = (exp(2*gammahi)-1)./(exp(2*gammahi)+1) - (((exp(2*gammahi)).^2-1).*yt2 - 2*ytp.*((exp(2*gammahi)).^2+1))./(4*vary'.*exp(2*gammahi))-h1;
     d2 = 4*exp(2*gammahi)./((exp(2*gammahi)+1).^2) - ((exp(2*gammahi)).^2 .*(yt2-2*ytp) + yt2 + ytp)./(2*vary'.*exp(2*gammahi)) - h2;
      ngammahi = gammahi - d1./d2;
      if sum(isnan(ngammahi))>0
          break
      end
      xidiff2 = sum(abs(ngammahi-gammahi));
      gammahi = ngammahi;
       end
       gamma(i,j,t,:)=gammahi;              
     Q = xi*xi';
     Q(isnan(Q)) = 0;
    sub_con(:,:,t) = sub_con(:,:,t) +Q;   %sub_con used for clustering
        end
    for tt=1:T
        subval = squeeze(gamma(i,j,tt,:));
        Hval = squeeze(gammaH(i,j,tt,:));
    [~,subind]=min(abs(subval-Hval')');
    subgroup(i,j,tt,:)=subind;
    end
end 
end  

% Update GammaH
for h=1:5
    for i=1:(V-1)
        for j=(i+1):V
            for t=1:T
whjlt(i,j,t,h) =  sqrt(sum(xiH(i,j,t,:,h))/(2*siggamma));
gammabar(i,j,t,h) = whjlt(i,j,t,h)*sum(gamma(i,j,t,:).*xiH(i,j,t,:,h))/(sum(xiH(i,j,t,:,h)));
            end
        end
    end
end
gammabar(isnan(gammabar)) = 0;
Mx = zeros(T,T);

for h=1:H
for i=1:(V-1)
    for j=(i+1):V
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
        if sum(abs(gamma2)<1e-3) >= (T-60)
        jum = [];
        jum2 = [1,jum'];
        jumval = 0;
        ss=1;
        res = 0;
        else
            [B] = glmnet(Mx,gamma2); %glmnet is used here
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
    end
gammaH(i,j,:,h)=res; 
    end
end    
gammaH(gammaH<1e-3) = 0;    
end

%update of siggamma
asigga = 9+0.25*squeeze((sum(sum(sum(sum(xiH))))));
siggamma = zeros(H,1);
for h=1:H
bsigga = 1+0.25*squeeze((sum(sum(sum(nansum(xiH(:,:,:,:,h).*(gamma-gammaH(:,:,:,h)).^2))))));
siggamma(i) = gamrnd(asigga(h),bsigga);
end
siggamma = 0.1;
%target function
for sub=1:nsub
    for t=1:T
   Yi = Yall(t,:,sub);
   yitjt(:,:,t,sub) = Yi'*Yi;
   yit2(:,:,t,sub)=Yi.^2 + Yi.^2';
    end
end
obj1 = -0.5*log(1-((exp(2*gamma)-1)./(exp(2*gamma)+1)).^2)-(yit2 - 2*((exp(2*gamma)-1)./(exp(2*gamma)+1)).*yitjt)./(2*varysub.*(1-((exp(2*gamma)-1)./(exp(2*gamma)+1)).^2));
obj1(obj1==-Inf) = 0;
obj2 = sum(sum(sum(nansum(obj1))));
obj4 = 0;
for h=1:H
   xi2 = xiH(:,:,:,:,h);
   psi2 = psiH(:,:,:,:,h);
   obj3 = -0.5*(xi2.*(gamma-gammaH(:,:,:,h)).^2)/(siggamma)- 0.5*gamma*log(siggamma)+ xi2.* log(psi2);% - xi2.*log(1-psi2);
   obj3(obj3==NaN) = 0;
   obj3(abs(obj3)==Inf) = 0;
   obj4 = obj4 + sum(sum(sum(nansum(obj3))));
end
obj5 = 0;
lambda = 1;
obj5 =obj5 - lambda*sum(sum(sum(nansum(gammaH(:,:,2:T,:)-gammaH(:,:,1:(T-1),:))))) - (a1+1)*log(siggamma)-b1/siggamma ;
objfun(niters) = obj2 + obj4 + obj5;

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

dy_conn = zeros(V,V,T,nsub);
for sub=1:nsub
   for t=1:T
       groupid = (sub_kmeans(:,t) == sub_kmeans(sub,t));
       dy_conn(:,:,t,sub) = mean(gamma(:,:,t,groupid),4);
   end
end
dy_conn(isnan(dy_conn)) = 0;

Graph = dy_conn;
%% Network CP
dy_cp = zeros(nsub,T);
dy_cp2 = zeros(V,V,T,nsub);
for sub=1:nsub
    [sub]
   for i=1:(V-1)
       for j=(i+1):V
        ploti = squeeze(dy_conn(i,j,:,sub))';    %change the input T data
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
maxdisti = 2;
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

