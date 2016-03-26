%%
%%Creat new network with dummy links
global incidenceFull;
global EstimatedTime;
global TurnAngles;
global LeftTurn;
global Uturn;   
global Scale; 

global incidenceFull_DM;
global EstimatedTime_DM;
global TurnAngles_DM;
global LeftTurn_DM;
global Uturn_DM; 
global ConstantAtts;
global Obs;
global prevLink;
global nDlinks; % number of dummy links
global nRlinks; % number of real links

globalVar;
%clear all;
%load('WSRL.mat');
%importfile('../LSatt.mat');


%% Counting number of dummy links
[nRlinks,m] = size(incidenceFull);
nDlinks = 0;
for k = 1: nRlinks
    nDlinks = nDlinks + size(find(incidenceFull(k,:)),2);
end
nDlinks

incidenceFull_DM = incidenceFull;
incidenceFull_DM(nRlinks+nDlinks, m + nDlinks) = 0;
EstimatedTime_DM = EstimatedTime;
EstimatedTime_DM(nRlinks+nDlinks, m + nDlinks) = 0;
TurnAngles_DM = TurnAngles;
TurnAngles_DM(nRlinks+nDlinks, m + nDlinks) = 0;
LeftTurn_DM = LeftTurn;
LeftTurn_DM(nRlinks+nDlinks, m + nDlinks) = 0;
Uturn_DM = Uturn;    
Uturn_DM(nRlinks+nDlinks, m + nDlinks) = 0;
prevLink = zeros(nRlinks+nDlinks,1);

for i= nRlinks+nDlinks +1: m + nDlinks
    incidenceFull_DM(:,i) = incidenceFull_DM(:,i - nDlinks);
    incidenceFull_DM(:,i - nDlinks) = 0;
end
ConstantAtts = incidenceFull_DM;
SavedIncidence = incidenceFull_DM;
% Load Scale
getScale();
%Scale_DM = Scale;
Scale(nRlinks+nDlinks + 1,size(Scale,2)) = 0;
startL = nRlinks;
for k = 1: nRlinks
    nextL = find(incidenceFull_DM(k,:));
    nd =  size(nextL,2); % number of dummies to add
    %Scale(k,:) = 0;
    for j = 1:nd
        a = nextL(j);
        incidenceFull_DM(k,a) = 0;
        ConstantAtts(k,a) = 0;
        incidenceFull_DM(k,startL+j) = 1;
        ConstantAtts(k,startL+j) = 0;
        prevLink(startL+j) = k;
        if (a <= nRlinks+nDlinks + 1)
            Scale(startL+j,:) = Scale(a,:);
        end
        for t = 1:nd
            incidenceFull_DM(startL+j, nextL(t)) = 1;
            ConstantAtts(startL+j, nextL(t)) = log(1/nd);
            EstimatedTime_DM(startL+j, nextL(t)) = EstimatedTime_DM(k,nextL(t)); 
            TurnAngles_DM(startL+j, nextL(t)) = TurnAngles_DM(k,nextL(t)); 
            LeftTurn_DM(startL+j, nextL(t)) = LeftTurn_DM(k,nextL(t)); 
            Uturn_DM(startL+j, nextL(t)) = Uturn_DM(k,nextL(t)); 
        end
    end
    for j = 1:nd
        a = nextL(j);
        EstimatedTime_DM(k,a) = 0; 
        TurnAngles_DM(k,a) = 0; 
        LeftTurn_DM(k,a) = 0; 
        Uturn_DM(k,a) = 0;
    end
    startL = startL + nd;
end

%% Mapping Observations
global nbobs;
maxLengthObs = size(Obs,2);
for iObs = 1:nbobs
    maxLengthObs = size(find(Obs(iObs,:)),2);
    Obs(iObs,1) = Obs(iObs,1) + nDlinks;
    Obs(iObs,maxLengthObs) = Obs(iObs,maxLengthObs) + nDlinks;
end

%% Transfer the KS attribute
global LSatt_DM;
global LSatt;
global isLinkSizeInclusive;
global Op;
if isLinkSizeInclusive ==  true
    LSatt_DM = objArray(nbobs);
    for iObs = 1:nbobs
        iObs
        startL = nRlinks;
        LSatt_DM(iObs).value =  LSatt(iObs).value;
        LSatt_DM(iObs).value(nRlinks + nDlinks, m + nDlinks) = 0;
        for k = 1: nRlinks
            nextL = find(SavedIncidence(k,:));
            nd =  size(nextL,2); % number of dummies to add
            for j = 1:nd
                a = nextL(j);
                for t = 1:nd
                    LSatt_DM(iObs).value(startL+j, nextL(t)) = LSatt_DM(iObs).value(k,nextL(t));
                end
                %LSatt_DM(iObs).value(k,a) = 0; 
            end
            for j = 1:nd
                a = nextL(j);
                LSatt_DM(iObs).value(k,a) = 0; 
            end
            startL = startL + nd;
        end
    end
end