%%   Example.m  %%%%%%%%%

% Description
% This example.m provide some examples about how to implement the functions
% 
% Fit a Bayesian product mixture model that imposes independent mixture priors 
% at each time scan, which results in time-varying clusters of samples designed to pool information.

% BPMM-PC.m   Bayesian Product Mxiture Modeling without covariates for pairwise correlation
% BPMM-PC.m   Bayesian Product Mxiture Modeling without covariates for precision matrix
% idPAC.m     Bayesian Product Mxiture Modeling with covariates for pairwise correlation
% idPMAC.m    Bayesian Product Mxiture Modeling with covariates for precision matrix

% Sample Data set
% example.mat: Nsub = 40 subjects are divided into 4 groups (1-10, 11-20, 21-30, 31-40).T = 300 scans, with V = 20
%              ROIs. The first 20 subjects have 3 changes points, the last
%              20 subjects have 4 change points. There are 2 binary
%              covariates divided the subjects into 4 groups
%              Y          fMRI scanning data    T*V*Nsub
%              X          covariates of all subjects   Nsub*2
%              changeT    True change points of each subjects Nsub*4
%              OmegaTall  True precision matrix of each group of subjects
%                         V*V*6*4
%              SigTrueall True pariwise correlation matrix of each group of
%                         subjects V*V*6*4 
% Output: 
% Graph         the estiamted pairwise correlation (V*V*T*Nsub)
% CP_all        change points for all subjects (Nsub*20)
% Clusteri      estimated cluster of all subjects 
% clusterCP     change points for all cluster (Nc*20)

load(example.mat)
%% BPMM-PC.m
% Assume we know that there are 4 groups of subjects, with H=4, and
% ncluster =4
[Graph,CP_all,Clusteri] = BPMM_PC(Y,4,4);

% Assume we do not have the prior information about the number of clusters
[Graph,CP_all,Clusteri] = BPMM_PC(Y);

% Once we get the cluster information, we could have a better estimation CP
[clusterCP] = clusterCP(Graph,Clusteri);

%% BPMM-PR.m

% Assume we know that there are 4 groups of subjects, with H=4, and
% ncluster =4
[Graph,CP_all,Clusteri] = BPMM_PR(Y,4,4);

% Assume we do not have the prior information about the number of clusters
[Graph,CP_all,Clusteri] = BPMM_PR(Y);

% Once we get the cluster information, we could have a better estimation CP
[clusterCP] = clusterCP(Graph,Clusteri);



%% idPAC.m
% Assume we have addition covaraites information
% Assume we know that there are 4 groups of subjects, with H=4, and
% ncluster =4
[Graph,CP_all,Clusteri] = idPAC(Y,X,H,ncluster);

% Assume we do not have the prior information about the number of clusters
[Graph,CP_all,Clusteri] = BPMM_PR(Y);

% Once we get the cluster information, we could have a better estimation CP
[clusterCP] = clusterCP(Graph,Clusteri);


%% idPMAC.m
% Assume we have addition covaraites information
% Assume we know that there are 4 groups of subjects, with H=4, and
% ncluster =4
[Graph,CP_all,Clusteri] = idPMAC(Y,X,H,ncluster);

% Assume we do not have the prior information about the number of clusters
[Graph,CP_all,Clusteri] = BPMM_PR(Y);

% Once we get the cluster information, we could have a better estimation CP
[clusterCP] = clusterCP(Graph,Clusteri);



