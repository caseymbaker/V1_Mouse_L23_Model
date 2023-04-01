%% code to select random pairs of neurons from the same ensemble

w = load('Gfile.mat').A; % load digraph from MakeNetwork.m 
w(w<3) = 0; %only keep strong connections

ne = 3200; % number of excitatory neurons
G = digraph(w(1:ne,1:ne)); %graph of only excitatory neurons

con = zeros(ne,ne); %initialize binary matrix of connections
%con is a matrix where each row is an excitatory neurons and
%the value is 1 if that neuron i is connected to neuron j in either
%direction
for i=1:ne
j = predecessors(G,i);
k = successors(G,i);
con(i,j) = 1; %mark presynaptic connections
con(i,k) = 1; %mark postsynaptic connections
end

%% cluster neurons into densely connected "ensembles"
mean_wc_dist = zeros(1,100); %vector of average within-cluster sums of point-to-centroid distances for each k
med_clust_sze = zeros(1,100); %vector of median cluster size for each k
for i= 1:100 
[idx,cat,sumd] = kmeans(con,i);
mean_wc_dist(i) = mean(sumd); 
sze = [];
for j = 1:max(idx)
    a = sum(idx == j);
    sze = [sze a];
end
    med_clust_sze(i) = median(sze);
end

% determined that k = 50 gives ensemble size ~ 30
%%
[idx,cat,sumd] = kmeans(con,50); %ensemble size should be ~ 30 neurons
bc = histogram(idx,length(unique(idx))).BinCounts;
histogram(bc(bc<900)) %histogram of ensemble sizes (excluding clusters > 900 neurons)
%find clusters of size <50 
clust = find(bc<50);
%%
c = randsample(clust,1);%picking a cluster at random
nrnidxs = find(idx==c);%get neurons from that cluster
save('clusterinfo.mat','c','idx','clust')

%% increase weights between cluster neurons
load('ExcData_old.mat')
for i=1:length(nrnidxs)
n1 = find(dataE(:,1) == nrnidxs(i));
n2 = find(ismember(dataE(:,2), nrnidxs));
n3 = intersect(n1,n2);
dataE(n3,3) = dataE(n3,3)*3.5;
end
save('ExcData.mat','dataE')
%% run if picking 5 sets of 2 neurons at random from the same cluster 
for i = 1:5
nodes = randsample(nrnidxs,2);
txt = ['nods' num2str(i) '.mat'];
save(txt,'nodes')
end
%% run if picking 5 sets of 2 random neurons
a = 0:3199;
for i = 1:5
nodes = randsample(nrnidxs,2);
txt = ['nods' num2str(i) '.mat'];
save(txt,'nodes')
end
