function [PCC, Latency, NetworkParams] = AnalyzeNetworkActivity(numneuron, numruns, folder)

%% inputs
% numneuron = number of neuron pairs stimulated (demo is 5)
% numruns = number of runs per pair (demo is 5)
% folder = folder with clusterinfo.mat, nods.mat, spikes.mat, and volt.mat,
% Gfile.mat
%% outputs 
% PCC = (1 x numneuron) vector with the PCC value for each neuron pair 

% Latency = (clustsize * numneuron) cell that gives the latency values for
% each neuron for each trial

% NetworkParams = (numneuron * 4) matrix where each column is a network
% parameter for the given neuron pair
% column 1 = weighted degree
% column 2 = closeness centrality
% column 4 = intersection
% column 5 = union
%% initializing variables 
numruns = numruns * 8; % *8 because there are 8 stims in the 250ms period
pertrial_meanV = zeros(numneuron,numruns); %mean voltage in 10ms before stim for each ensemble neuron for each trial
pertrial_medV = zeros(numneuron,numruns); %median voltage in 10ms before stim for each ensemble neuron for each trial
ensemble = zeros(numneuron,numruns); %number of active neurons in 10ms after stim for each trial
pertrialVoltage = cell(numneuron,numruns); %cell of voltage traces for each ensemble neuron for each trial
spks = cell(numneuron,numruns); %cell of spikes for each ensemble neuron for each trial

%cd(folder)
load('clusterinfo.mat')
cluster = find(idx == c);
sz = length(cluster); %number of neurons in ensemble
%% calculating trial info
for i = 1:numneuron
    i
z = 1;

load(['spikes' num2str(i) '.mat']);
load(['volt' num2str(i) '.mat']);
for k = 1:numruns/8

    
    v = squeeze((volts(k,cluster,:)));
    s = (k-1)*4000+1;
    a = arr(s:(s+3199),:);
    s = a(cluster,:);
    start = 10000;
    for n = 1:8
        spi = s(:,(start-1):(start+100));
        vols = v(:,start-100:start-1);
        medv = median(vols(:));
        meanv = mean(vols(:));
        pertrial_meanV(i,z) = meanv;
        pertrial_medV(i,z) = medv;
        ensemble(i,z) = sum(sum(spi)>0);
        pertrialVoltage{i,z} = v(:,start-100:start+100);
        spks{i,z} = spi;
        start = start + 333;
        z=z+1;
    end
end
end

%% Calculating PCC

n = 5; % n*3 = number of bins 
range = linspace(min(pertrial_meanV(:)),max(pertrial_meanV(:)),n);
range2 = linspace(min(pertrial_meanV(:)),max(pertrial_meanV(:)),n)+(range(2)-range(1))/3;
range3 = linspace(min(pertrial_meanV(:)),max(pertrial_meanV(:)),n)+2*(range(2)-range(1))/3;

ratio = zeros(numneuron,(n-1)*3);

for i = 1:numneuron
    thistrial = ensemble(i,:);
    yestrials = find(thistrial>=0.75*sz);
    notrials = find(thistrial< 0.75*sz);
    bcy = histogram(pertrial_meanV(i,yestrials),range).BinCounts;
    bcn = histogram(pertrial_meanV(i,:),range).BinCounts;
    ratio(i,1:3:(n-1)*3-2) = bcy./bcn;
    bcy = histogram(pertrial_meanV(i,yestrials),range2).BinCounts;
    bcn = histogram(pertrial_meanV(i,:),range2).BinCounts;
    ratio(i,2:3:(n-1)*3-1) = bcy./bcn;
    bcy = histogram(pertrial_meanV(i,yestrials),range3).BinCounts;
    bcn = histogram(pertrial_meanV(i,:),range3).BinCounts;
    ratio(i,3:3:(n-1)*3) = bcy./bcn;
end

%voltages
x = zeros(1,(n-1)*3);
x(1,1:3:(n-1)*3-2) = range(1:n-1);
x(1,2:3:(n-1)*3-1) = range2(1:n-1);
x(1,3:3:(n-1)*3) = range3(1:n-1);

ensemblelog = logical(ensemble>0.75*sz);
activ = sum(ensemblelog,2)/size(ensemble,2);
col = (activ-min(activ))/(max(activ)-min(activ));

auc = zeros(1,numneuron);
for i = 1:numneuron

plot(x(1:8),ratio(i,1:8),'Color',[(1-col(i)) 0 col(i)],'LineWidth',2)
xlabel('Pre-Stim Voltage')
ylabel('Probability of Activating Ensemble')
xlim([min(range),max(range)])
hold on;
auc(i) =  trapz(x(1:8),ratio(i,1:8));

end

r2 = zeros(1,numneuron);
vpt5 = zeros(1,numneuron);
for i=1:numneuron
x1 = x(1:8);
y1 = ratio(i,1:8);
ft = fittype( '1/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1 1];
% Fit model to data.
[fitresult, gof] = fit( x1', y1', ft, opts );
xnew = linspace(-0.08,-0.055,100);
y = fitresult(xnew);
r2(i) = gof.rsquare;
if (max(y)>0.05)
vpt5(i) = xnew(find(y>0.05,1));
else
vpt5(i) = NaN;
end
s = scatter(x1,y1,'filled')
hold on;
plot(fitresult)
end

PCC = vpt5;

%% Get latency measure

latency = cell(sz,numneuron);
f = find(idx == c);
for k = 1:numneuron
y = find(ensemble(k,:)>25);
for j = 1:length(y)
    d = y(j);
    v = spks{k,d};
    for i = 1:sz
        latency{i,k} = [latency{i,k} find(v(i,:),1)];
    end        
end
end

%remove latency values for stimulated neurons
for i = 1:numneuron
  txt =  ['nods' num2str(i) '.mat'];
  load(txt)
  a = find(f == nodes(1));
  b = find(f == nodes(2));
  latency{a,i} = NaN;
  latency{b,i} = NaN;
end


%% get network params
clustnrns = find(idx == c);
load('Gfile.mat')
newG = A(clustnrns,clustnrns);
newG(newG<3) = 0;
[X,Y] = meshgrid(1:sz,1:sz);
Y = Y(:);
X = X(:);
weighted_deg = sum(newG,2);
newG = newG(:);
pos = find(newG > 0);
wt = newG(pos);
cost = max(wt)-wt+0.01;
G = digraph(Y(pos),X(pos),newG(pos));
presyn = cell(2,numneuron);
postsyn = cell(2,numneuron);
outclose = centrality(G,'outcloseness','Cost',cost);
outdeg_weighted = zeros(2,numneuron);
centrality_outcloseness = zeros(2,numneuron);
samesuccess = zeros(1,numneuron);
numneuronspost = zeros(1,numneuron);
for i = 1:numneuron
  txt =  ['nods' num2str(i) '.mat'];
  load(txt)
  a = find(clustnrns == nodes(1));
  b = find(clustnrns == nodes(2));  
  samesuccess(i) = length(intersect(successors(G,a), successors(G,b)));
  numneuronspost(i) = length(union(successors(G,a), successors(G,b)));
  outdeg_weighted(1,i) = weighted_deg(a);
  outdeg_weighted(2,i) = weighted_deg(b);
  centrality_outcloseness(1,i) = outclose(a);
  centrality_outcloseness(2,i) = outclose(b);

end

NetworkParams = [sum(outdeg_weighted)',sum(centrality_outcloseness)',samesuccess', numneuronspost'];
NetworkParams = normalize(NetworkParams,1);


end
