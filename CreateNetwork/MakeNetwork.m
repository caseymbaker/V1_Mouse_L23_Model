%% Number of Neurons in each type
% 80% excitatory, 20% inhibtory 
NE = 3200;   % number of excitatory neurons
NI = 800;   % number of inhibitory neurons
N = NE+NI;  % total number of neurons

% breakdown of inhibitory neurons
NPV = 330; % number of PV neurons
NSOM = 330; % number of SOM neurons
NVIP = 140; % number of VIP neurons

%% Background Current (used at the end of code)
bc = 75 % if doing stimulation/increasing ensemble connection weights
bc = 100 % if looking at spontaneous activity

numframs = 20000; % if doing stimulation
numframs = 50000; % if doing spontaneous activity
%% Preferred Orientation (po)
po_exc = linspace(0,pi,NE)+ 1/10*(rand(1,NE)-.5); %goes from 0 to pi with some variability 
save('po_exc.mat','po_exc')

po_PV = linspace(0,pi,NPV); % po PV 
po_SOM = linspace(0,pi,NSOM); % po SOM
po_VIP = linspace(0,pi,NVIP); % po VIP 
po_all = [po_exc, po_PV, po_SOM, po_VIP]; 

ab = pdist(po_all'); % po distance between each neuron (smaller number = more similar)
ab(ab>pi/2) = pi-ab(ab>pi/2); % correcting dist measure so that perpendicular measures are max distance
ab = squareform(ab)+ 1/20*(rand(N,N)-.5); %adding random noise to similarity metric
ab = pi-ab; %convert distance to similarity (larger number = more similar)
ab = ab - diag(diag(ab)); %remove self-connections

%% Locations 
% all neurons are randomly distributed 
% mouse RF similarity is not related to cortical distance (Cossell et al.,2015, Nature)
x = rand([1,N]);
y = rand([1,N]);
loc = [x;y]; %neuron locations in a 2-D plane
dis = squareform(pdist(loc')); %euclidean distance between neurons
m = max(dis(:))+1; 
dis = (m-dis)/m; %normalizing the distances
dis = dis - diag(diag(dis)); %remove self-connections

%% Exc-Exc Connectivity 
% all connections have a lognormal distribution
wEE = ab(1:3200,1:3200).*(.25*dis(1:3200,1:3200)); %weight bias 
%weight of connection is based mostly on po similarity and a little on cortical distance
%rows = presynaptic neuron, columns = postsynaptic neuron
wEE2 = ab(1:3200,1:3200); %initializing final weight matrix 

for i = 1:3200 %loop through each presynaptic neuron
    aa = wEE(:,i);
    tokeep1 = prctile(aa(:),85); %keep 15% of connections
    v = size(aa(aa>tokeep1),1); %number of postsynaptic neurons for neuron i
    rnum = lognrnd(-0.25,0.75,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order lognormal weights
    [out2 ,idx2] = sort(aa,'descend'); %order weight biases
    c = idx2(1:v); 
    d = idx2(v+1:length(idx2));
    wEE2(c,i) = out; %assign connection weights based on weight bias
    wEE2(d,i) = 0; %85th percentile of connections and below = no connection
end

bins = 10.^(-1:.2:1.2);
wEE2 = 0.3021*wEE2; %converting weight (unitless) to EPSP (mV) (based on weight/EPSP relationship found using Brian2.0)
wEE2(wEE2>10) = 10; %threshold EPSP
h = histogram(wEE2(wEE2>0),bins) %visualize weight distribution
set(gca,'xscale','log')
xlabel('EPSP Amplitude (mV)')
ylabel('Number of Connections')
title('E->E')
h.FaceColor = [0 0 0];
EEv = wEE2; %mV EPSP
wEE = wEE2/0.3201; %weight
%% E->PV connectivity (repeat above steps for rest of connection types)

wEPV = ab(1:3200,3201:3530);
wEI2 = ab(1:3200,3201:3530);

for i = 1:3200
    aa = wEPV(i,:);
    tokeep1 = prctile(aa(:),30); %keep 70% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(1.25,0.85,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); 
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wEI2(i,c) = out';
    wEI2(i,d) = 0;
end
bins = 10.^(-1:.2:1.2);
wEI2 = 0.3021*wEI2;
wEI2(wEI2>10) = 10;
max(wEI2(:))
h = histogram(wEI2(wEI2>0),bins)
set(gca,'xscale','log')
xlabel('EPSP Amplitude (mV)')
ylabel('Number of Connections')
title('E->PV')
h.FaceColor = [0 0 0];
wEPV = wEI2/.3021; %weight (unitless)
EPVv = wEI2; %EPSP (mV)
%% PV->E  connectivity

%PV neurons synapse strongly onto EXC neurons that synapse strongly onto
%them
wIE = ab(3201:3530,1:NE);
wIE2 = ab(3201:3530,1:NE);

for i = 1:330
    aa = wIE(i,:);
    tokeep1 = prctile(aa(:),20);%keep 80% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(-0.25,0.85,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); 
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wIE2(i,c) = out';
    wIE2(i,d) = 0;
end
bins = 10.^(-1:.1:1.2);
wIE2(wIE2>10) = 10;
max(wIE2(:))
h = histogram(wIE2(wIE2>0),bins)
set(gca,'xscale','log')
xlabel('IPSP Amplitude (mV)')
ylabel('Number of Connections')
title('PV -> E')
h.FaceColor = [0 0 0];
wPVE = wIE2/.6537; %weight %converting weight (unitless) to EPSP (mV) (based on weight/IPSP relationship found using Brian2.0)
PVEv = wIE2; %IPSP
%% E --> SOM 

wESOM = ab(1:NE,3531:3860); 
wEI2 = ab(1:NE,3531:3860); 

for i = 1:3200
    aa = wESOM(i,:);
    tokeep1 = prctile(aa(:),30); %keep 70% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(1.2,0.75,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); %
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wEI2(i,c) = out';
    wEI2(i,d) = 0;
end
bins = 10.^(-1:.2:1.2);
wEI2 = 0.3021*wEI2;
wEI2(wEI2>15) = 15;
max(wEI2(:))
h = histogram(wEI2(wEI2>0),bins)
set(gca,'xscale','log')
xlabel('EPSP Amplitude (mV)')
ylabel('Number of Connections')
title('E->SOM')
h.FaceColor = [0 0 0];
wESOM = wEI2/.3021; %weight, histogram above is in mV
ESOMv = wEI2; %EPSP
%% SOM-->E

wSOME = dis(3531:3860,1:3200); %connection based on cortical distance 
wIE2 = dis(3531:3860,1:3200);

for i = 1:330
    aa = wSOME(i,:);
    tokeep1 = prctile(aa(:),20);%keep 80% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(0,0.25,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); %
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wIE2(i,c) = out';
    wIE2(i,d) = 0;
end
bins = 10.^(-1:.1:1.2);
wIE2(wIE2>10) = 10;
wIE2 = wIE2/2;
h = histogram(wIE2(wIE2>0),bins)
set(gca,'xscale','log')
xlabel('IPSP Amplitude (mV)')
ylabel('Number of Connections')
title('SOM->E')
h.FaceColor = [0 0 0];
wSOME = wIE2/.6537; %weight
SOMEv = wIE2; %IPSP
%% PV --> PV 
wPVPV = dis(3201:3530,3201:3530); 
wPVPV2 = dis(3201:3530,3201:3530); 

for i = 1:330
    aa = wPVPV(:,i);
    tokeep1 = prctile(aa(:),20);%keep 80% of connections
    v = size(aa(aa>tokeep1),1);
    rnum = lognrnd(-.4,0.3,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); %
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wPVPV2(c,i) = out;
    wPVPV2(d,i) = 0;
end
bins = 10.^(-1:.1:1.2);
wPVPV2(wPVPV2>10) = 10;
max(wPVPV2(:))
h = histogram(wPVPV2(wPVPV2>0),bins)
set(gca,'xscale','log')
xlabel('IPSP Amplitude (mV)')
ylabel('Number of Connections')
title('PV-PV')
h.FaceColor = [0 0 0];
wPVPV = wPVPV2/.6537; %weight
PVPVv = wPVPV2; %IPSP
%% SOM --> PV 

wSOMPV = dis(3531:3860,3201:3530); 
wIE2 = dis(3531:3860,3201:3530);

for i = 1:330
    aa = wSOMPV(i,:);
    tokeep1 = prctile(aa(:),30);%keep 80% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(-.4,0.3,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); %
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wIE2(i,c) = out';
    wIE2(i,d) = 0;
end
bins = 10.^(-1:.1:1.2);
wIE2(wIE2>10) = 10;
wIE2 = wIE2/3;
h = histogram(wIE2(wIE2>0),bins)
set(gca,'xscale','log')
xlabel('IPSP Amplitude (mV)')
ylabel('Number of Connections')
title('SOM->PV')
h.FaceColor = [0 0 0];
wSOMPV = wIE2/.6537; %weight
SOMPVv = wIE2; %IPSP
%% E --> VIP

wEVIP = ab(1:NE,3861:4000); 
wEI2 = ab(1:NE,3861:4000); 

for i = 1:3200
    aa = wEVIP(i,:);
    tokeep1 = prctile(aa(:),30); %keep 70% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(1,.5,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); %
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wEI2(i,c) = out';
    wEI2(i,d) = 0;
end
bins = 10.^(-1:.2:1.2);
wEI2 = 0.3021*wEI2;
wEI2(wEI2>15) = 15;
max(wEI2(:))
h = histogram(wEI2(wEI2>0),bins)
set(gca,'xscale','log')
xlabel('EPSP Amplitude (mV)')
ylabel('Number of Connections')
title('E->VIP')
h.FaceColor = [0 0 0];
wEVIP = wEI2/.3021; %weight, histogram above is in mV
EVIPv = wEI2; %EPSP
%% VIP --> SOM

wVIPSOM = ab(3861:4000,3531:3860); %rand(3200,330);
wEI2 = ab(3861:4000,3531:3860); %rand(3200,330);

for i = 1:140
    aa = wVIPSOM(i,:);
    tokeep1 = prctile(aa(:),40); %keep 60% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(-.4,0.3,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); %
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wEI2(i,c) = out';
    wEI2(i,d) = 0;
end
bins = 10.^(-1:.2:1.2);
wEI2 = 0.6537*wEI2;
wEI2(wEI2>15) = 10;
wEI2 = wEI2;
max(wEI2(:))
h = histogram(wEI2(wEI2>0),bins)
set(gca,'xscale','log')
xlabel('IPSP Amplitude (mV)')
ylabel('Number of Connections')
title('VIP -> SOM')
h.FaceColor = [0 0 0];
wVIPSOM = wEI2/.6537; %weight, histogram above is in mV
VIPSOMv = wEI2; %EPSP
%% SOM --> VIP

wSOMVIP = dis(3531:3860,3861:4000);
wIE2 = dis(3531:3860,3861:4000);

for i = 1:330
    aa = wSOMVIP(i,:);
    tokeep1 = prctile(aa(:),20);%keep 80% of connections
    v = size(aa(aa>tokeep1),2);
    rnum = lognrnd(0,0.25,v,1); %generate lognorm distribution
    [out ,idx] = sort(rnum,'descend'); %order of rand log dist
    [out2 ,idx2] = sort(aa,'descend'); %
    c = idx2(1:v);
    d = idx2(v+1:length(idx2));
    wIE2(i,c) = out';
    wIE2(i,d) = 0;
end
bins = 10.^(-1:.1:1.2);
wIE2(wIE2>10) = 10;
wIE2 = wIE2/2;
h = histogram(wIE2(wIE2>0),bins)
set(gca,'xscale','log')
xlabel('IPSP Amplitude (mV)')
ylabel('Number of Connections')
title('SOM->VIP')
h.FaceColor = [0 0 0];
wSOMVIP = wIE2/.6537; %weight
SOMVIPv = wIE2; %IPSP
%% Create matrix of weights between all 4000 neurons

%not connected 
wPVSOM = zeros(NPV,NSOM);
wPVVIP = zeros(NPV,NVIP);
wSOMSOM = zeros(NSOM,NSOM);
wVIPE = zeros(NVIP,NE);
wVIPVIP = zeros(NVIP,NVIP);
wVIPPV = zeros(NVIP,NPV);

%weights based on each presynaptic neuron type
wE = [wEE,wEPV,wESOM,wEVIP];
wPV = [wPVE,wPVPV,wPVSOM,wPVVIP];
wSOM = [wSOME,wSOMPV,wSOMSOM,wSOMVIP];
wVIP = [wVIPE,wVIPPV,wVIPSOM,wVIPVIP];

w = [wE;wPV;wSOM;wVIP]; % all weights

vE = [EEv,EPVv,ESOMv,EVIPv];
vPV = [PVEv,PVPVv,wPVSOM,wPVVIP];
vSOM = [SOMEv,SOMPVv,wSOMSOM,SOMVIPv];
vVIP = [wVIPE,wVIPPV,VIPSOMv,wVIPVIP];
v = [vE;vPV;vSOM;vVIP]; %all E/I PSPs (mV) 

%% Creating directed graph 

A=w; %create a copy of weight matrix
G = digraph(w);

de = G.outdegree()+G.indegree(); %total degree of each neuron

%remove neurons with few connections 
rm = find(de<=2); 
g2 = rmnode(G,rm);

A(:,rm) = [];
A(rm,:) = [];


[row,col] = find(A~=0); %row = presynaptic, column = postsynaptic
B = A(:); %vector of weights
B(B==0) = []; %remove non-connections

%to import to python: 
%row 1: presynaptic neuron, row 2: postsynaptic neuron, row 3: weight of
%connection
data = [row,col,B]; 

save('connections.mat','data') %save connection weights
save('Gfile.mat','A') %save digraph

% save excitatory and inhibitory connection info 
dataE = data(data(:,1,:)<=NE,:); %excitatory connections
dataI = data(data(:,1,:)>NE,:); %inhibitory connections 

save('ExcData.mat','dataE') %should be ExcData_old.mat if doing stimulation (see picknrns.m)
save('InhData.mat','dataI')
%% Create/Save Background Noise File

N = 4000;
NE = 3200;
memC = 200; %membrane capacitance
rstim = rand(N,numframs)*5*memC; %creating random noise between 0 and 5 mV and converting to current
rstim2 = randi([-1 1],N,numframs); %making some of the noise hyperpolarizing and some depolarizing
%making noise sparse 
rstim3 = randi([0 1],N,numframs);
for i = 1:5 
    rstim3 = rstim3.*randi([0 1],N,numframs);
end
rstim4 = rstim.*rstim2.*rstim3;

rstim4(1:NE,:) = rstim4(1:NE,:)+bc; %exc neurons have some background current to maintain their firing rate and resting potential
imagesc(rstim4)
save('stimfile.mat','rstim4')


