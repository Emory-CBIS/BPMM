function [CP_all] = clusterCP(Gall,cluster)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clusterCP.m: 
% BPMM (Bayesian Product Mxiture Modeling)
% Description: Estimating change points based on the cluster, given either
%              dynamic pairwise correlation or dynamic partial correlation
%              of subjects. This function is a follow-up function based on
%              the result of BPMM_PC.M, BPMM_PR.M, idPAC.M, or idPMAC.M.
%              Assume that all subjects within one culster share the same
%              number and location of change points.
%
% Usage:
% [CP_all] = clusterCP(Gall,cluster);
%               
%
% Input:
% Gall          The estimated pairwise correlation or partial correlation 
%               (V*V*T*Nsub), V is the number of ROIs, T is the number of
%               scans, Nsub is the number of subjects     
% cluster       estimated cluster of all subjects 
%
% Output: 
% CP_all        change points for all cluster (Nc*20), Nc is the number of
%               clusters defined in the input variable (cluster), 20 is the
%               max number of change points of each cluster
% 
% AUTHORS:
% Suprateek Kundu     Emory University
% Jin Ming            Emory University
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


V = size(Gall,1);
T = size(Gall,3);
nsub = size(Gall,4);

dy_cp = zeros(nsub,T);
dy_cp2 = zeros(V,V,T,nsub);

Nc = max(cluster);
for sub=1:nsub
   for i=1:(V-1)
       for j=(i+1):V
        ploti = squeeze(Gall(i,j,:,sub))';   
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

perceni = 0.4;
maxdisti = 1;
cp_est = zeros(Nc,20);
for sub=1:Nc
    cp_esti = sum(dy_cp(cluster==sub,:));
    cp_esti(T)=0;
    cp_estj = cp_clus(cp_esti,perceni,maxdisti);
    cp_est(sub,1:length(cp_estj)) = cp_estj;
end

end





