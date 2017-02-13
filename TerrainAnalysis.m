%% Terrain Analysis
% This script analyzes the surface geometry and creates a SWE model based
% on terrain parameters. Extrapolates SWE based on a stepwise linear regression.

%inputs: 
% dem geotiff UTM format
% extent geotiff format and snapped to dem

% this code assumes that the input files are in a projection
% that uses some sort of equal X,Y,Z units, like a UTM or albers coordinate system.
% 
% Louis Sass - written 2016.06.01
% last edited 2017.01.31

close all
clear all
tic
addpath data/  
addpath functions/
dbstop if error
opengl software;
w = waitbar(1/100, 'your computer is very slow :(');

%% set evaluation parameters here
radius = 370; % set radius over which to evaluate curvature, in the same units as the dem cellsize
distance = [10,100]; % set distance range for shelter term Sx, where minimum must be >= cellsize
near_distance = [10,20]; % set distance ranges for slopebreak term Sb
far_distance = [310,2000];
direction = 048; % set direction for prevailing wind direction for shelter term Sx 
dem = ('2015.08.13.WolvDEMreg_10m.tif');%DEM goes here
mask = ('2015mask_10m.tif');%Extent goes here
pick = ('wolverine_2015_Set11_constant.txt');%Pick file goes here

findOptimumCurvature = ('no'); % set to yes to find best curvature values, then turn off and re-run using those values above
findOptimumSx = ('no'); % set to yes to find best Sx values, then turn off and re-run using those values above
plots = ('off'); %turn on extra plots for troubleshooting
isHelicopter = ('yes'); %set to 'no' for ground data, 'yes' for mixed or helicopter data
maxslope = 40; %set the slope mask here to ditch steep terrain data

names = {'Z','Sb','Sx','curvature','northness','eastness','slope'}; %Possible names: 'Z', 'slope', 'aspect', 'northness', 'eastness', 'curvature','Sx','Sb'
labels = {'elevation [1000 m]' 'slopebreak' 'shelter' 'curvature' 'northness' 'eastness' 'slope [{\circ}]'};

%% import and format the data

[DEM.Z, DEM.ref] = geotiffread(dem); 
if DEM.ref.CellExtentInWorldX == DEM.ref.CellExtentInWorldY
    DEM.cellsize = DEM.ref.CellExtentInWorldX;
else
    error('DEM should be on a square grid')
end
DEM.Z(DEM.Z==DEM.Z(1,1))=NaN; %DEM.Z = the DEM value map
DEM.info = geotiffinfo(dem);
DEM.R = DEM.info.RefMatrix;
DEM.x = repmat((DEM.ref.XWorldLimits(1) + DEM.ref.CellExtentInWorldX/2:DEM.ref.CellExtentInWorldX:DEM.ref.XWorldLimits(2) - DEM.ref.CellExtentInWorldX/2),DEM.ref.RasterSize(1),1);
DEM.y = repmat(flipud((DEM.ref.YWorldLimits(1) + DEM.ref.CellExtentInWorldY/2:DEM.ref.CellExtentInWorldY:DEM.ref.YWorldLimits(2) - DEM.ref.CellExtentInWorldY/2)'),1,DEM.ref.RasterSize(2));

waitbar(1/10, w, 'DEM import complete, importing the mask');
[GL.on, GL.ref] = geotiffread(mask); %GL.on = glacier extent mask
if GL.ref.CellExtentInWorldX == GL.ref.CellExtentInWorldY
    GL.cellsize = GL.ref.CellExtentInWorldX;
else
    error('extent should be on a square grid')
end
GL.info = geotiffinfo(mask);
GL.R = GL.info.RefMatrix;

waitbar(2/10, w, 'importing GPR');
SWE.data = readtable(pick); %SWE.data = GPR pick data 
[SWE.data.x SWE.data.y] = projfwd(DEM.info,SWE.data.lat, SWE.data.long);
swemax = max(SWE.data.SWE); %Max pick value in GPR
waitbar(3/10, w, 'all imports complete');

%% check inputs

if GL.cellsize ~= DEM.cellsize 
    error('someone fucked up - inputs should have the same cellsize')
end  
        
%% evaluate terrain

waitbar(4/10, w, 'calculating terrain parameters');
[DEM.slope,DEM.aspect,DEM.eastness,DEM.northness] = CalcTerrainParams(DEM.Z,DEM.cellsize,w);
DEM.shade = 2.*DEM.eastness + DEM.northness;

%%
waitbar(5/10, w, 'calculating curvature');
[DEM.curvature] = Curvature(DEM.Z,DEM.cellsize,radius, 'method', 'xy'); %method 'xy' is generally much more useful.
waitbar(6/10,w,'calculating Sx')
temp_sx = Sx(DEM.Z, DEM.cellsize, distance(1),distance(2), direction);
%temp_sx(temp_sx>=0) = 0; %If Sx is supposed to be an erosional term then +values don't make sense. In practise it doesn't correlate well with snow depths if you restrict it to negative values 
[DEM.Sx] = temp_sx;
waitbar(7/10,w,'calculating Sb')
temp_sb = Sb(DEM.Z, DEM.cellsize, near_distance,far_distance, direction);
temp_sb(temp_sb<=0) = 0; % +slopebreaks are well correlated with drifts, -slopebreaks are not necessarily wind scoured, the correlation here is better if you restrict Sb to + values 
[DEM.Sb] = temp_sb;
waitbar(8/10,w, 'done with terrain, plotting....');

%% terrain parameter figures for the whole DEM

for n = 1:length(names)
figure ();
colormap(jet)
if strcmp(names{n},'curvature')==1
    clims = [-0.5 0.5];
    imagesc(DEM.(names{n}), 'alphadata', ~isnan(DEM.Z),clims);hold on
else
    imagesc(DEM.(names{n}), 'alphadata', ~isnan(DEM.Z));hold on
end
axis ij;
axis image;   
colorbar;
xlabel('east');
ylabel('north');
text(100, 100, labels{n});
end

close (w)

%% 
Xtemp = [DEM.x(:),DEM.y(:),DEM.Z(:)];
index = ~isnan(Xtemp(:,3));

if strcmp(isHelicopter,'yes')
    w = waitbar(3/100, 'correcting XY positions for Helicopter Data');
X = double(Xtemp(index,:));
T = delaunayn(X);
waitbar(5/100, w, 'your computer is pathetically slow');
Xi = [SWE.data.x,SWE.data.y,(SWE.data.elev + 12)];% if GPR output changes to ellipsoid heights then remove the +12
gprIndicies = dsearchn(X,T,Xi); 
else
     w = waitbar(3/10, 'extracting XY positions');
X = double(Xtemp(index,1:2));
T = delaunayn(X);
Xi = [SWE.data.x,SWE.data.y];
gprIndicies = dsearchn(X,T,Xi);   
end

waitbar(7/10, w, 'extracting terrain parameters');
for n = 1:length(names);
    waitbar(n/10, w, ['calculating ' (names{n}) ' at GPR data locations']);
    allvalues = DEM.(names{n});
    vectorvalues = allvalues(:);
    nonanvalues = vectorvalues(index,1);
    GPR.(names{n}) = nonanvalues(gprIndicies);
end

%% test plots

ind = find(GPR.slope<=maxslope);
SWE.c = ((SWE.data.x(ind) - DEM.ref.XWorldLimits(1))./DEM.cellsize) + 1/2; %interp
SWE.r = ((DEM.ref.YWorldLimits(2) - SWE.data.y(ind))./DEM.cellsize) + 1/2;

if strcmp(plots, 'on')==1
for n= 1:length(names)
figure ();
colormap(jet)
scatter(SWE.c,SWE.r,3,GPR.(names{n}))
axis ij;
axis image;
colorbar;
caxis([0,swemax]);
xlabel('east');
ylabel('north');
text(700, 900, labels{n});
end
end

%% regression
waitbar(8/10, w, 'regressing')
xySWE = SWE.data.SWE ./ cosd(GPR.slope); %corrects for slope
[NormalxySWE SWEmu SWEsigma] = zscore(xySWE); % standardized
TerrainTable = table(NormalxySWE(ind));
TerrainTable.Properties.VariableNames{'Var1'} = 'SWE';
for n= 1:length(names);
    [add mu sigma] = zscore(GPR.(names{n})); %standardized
    TerrainTable.(names{n}) = add(ind);
    muV.(names{n}) = mu;
    sigmaV.(names{n})= sigma;
end


%% Regression
%'Criterion' allows you to choose what type of cuttoff gets used for including params. choices are 'Variance', 'SSE','AIC','BIC', 'RSquared',or 'AdjRSquared'.
%'PEnter' allows you to choose the cutoff for the designated criterion
%'Upper' allows you to limit model complexity, where 'linear' means you will only include linear terms rather than interaction terms, etc.

mdl = stepwiselm(TerrainTable,'ResponseVar','SWE','Criterion','AdjRsquared','PEnter',0.01,'PRemove',-0.005,'Upper','linear','verbose',2); %,'PredictorVars',{'Z','Sb','curvature','eastness'},
%'Criterion','SSE','PEnter',0.05,'PRemove',0.10,
%'Criterion','AIC','PEnter',0,'PRemove',0.01,

%% Loop to find Curvature
if strcmp(findOptimumCurvature, 'yes')==1 
 wb = waitbar(1/100, 'OptimizeCurvature Loop');  
 [CurvatureValues] = OptimizeCurvature(DEM.Z, DEM.cellsize, TerrainTable, gprIndicies, ind, wb);  
end

%% Loop to find Sx
if strcmp(findOptimumSx, 'yes')==1 
  wb = waitbar(1/100, 'OptimizeSx Loop');  
 [SxValues] = OptimizeSx(DEM.Z, DEM.cellsize, TerrainTable, gprIndicies, ind, wb);   
end

%% Overlay residuals on hillshade

figure()
ax1 = axes;
imagesc(DEM.shade, 'alphadata', ~isnan(DEM.Z));hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
ax2 = axes;
linkaxes([ax1,ax2])
scatter(SWE.c,SWE.r,1,mdl.Residuals.Raw); hold on
axis ij;
axis image;
caxis([-swemax,swemax]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'jet')
set([ax1,ax2],'Position',[.10 .12 .7 .78]);
cb2 = colorbar(ax2,'Position',[.88 .12 .05 .78]);
title(ax1,'Model Residuals');
xlabel(ax1,'east [m]');
ylabel(ax1,'north [m]');

%% extrapolation

waitbar(9/10, w, 'extrapolating')
XSWE.alln = mdl.Coefficients.Estimate(1);
for n = 1:length(mdl.PredictorNames);
    predictor = mdl.PredictorNames{n};
XSWE.alln = XSWE.alln + mdl.Coefficients.Estimate(n+1) * ((DEM.(predictor)-muV.(predictor))./ sigmaV.(predictor));
end

XSWE.all = XSWE.alln .* SWEsigma + SWEmu; %returns it from standardized to actual
%% extrapolated SWE values

figure()
imagesc(XSWE.all, 'alphadata', ~isnan(DEM.Z));hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
h1 = gca;
caxis([0,swemax]);
colormap(h1,jet); hold on
colorbar('eastoutside');
xlabel('east');
ylabel('north');
text(700, 900, 'extrapolated SWE');

%% masking out the glacier data

[m,n] = size(GL.on);
GL.aligned = zeros(m,n);
GL.aligned(GL.on==0) = 1;
addx = (GL.ref.XWorldLimits(1) - DEM.ref.XWorldLimits(1))/DEM.ref.CellExtentInWorldX;
add = zeros(m,addx).*255;
GL.aligned = [add,GL.aligned];
addx = (DEM.ref.XWorldLimits(2) - GL.ref.XWorldLimits(2))/DEM.ref.CellExtentInWorldX;
add = zeros(m,addx).*255;
GL.aligned = [GL.aligned,add];
[m,n] = size(GL.aligned);
addy = (DEM.ref.YWorldLimits(2) - GL.ref.YWorldLimits(2))/DEM.ref.CellExtentInWorldY;
add = zeros(addy,n).*255;
GL.aligned = [add;GL.aligned];
addy = (GL.ref.YWorldLimits(1) - DEM.ref.YWorldLimits(1))/DEM.ref.CellExtentInWorldY;
add = zeros(addy,n).*255;
GL.aligned = [GL.aligned;add];
if size(GL.aligned)~=size(DEM.Z);
    error('the mask is not aligned with the data')
end
GL.area = nansum(nansum(GL.aligned));
XSWE.on = NaN(size(XSWE.all));
XSWE.on(GL.aligned==1) = XSWE.all(GL.aligned==1);
GL.SWE = nansum(nansum(XSWE.on))/GL.area;

%% Glacier SWE
figure()

imagesc(XSWE.all, 'alphadata', GL.aligned);hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
h1 = gca;
caxis([0,swemax]);
colormap(h1,jet); hold on
colorbar('eastoutside');
xlabel('east');
ylabel('north');
text(50, 50, 'extrapolated SWE');
text(50, 100, ['based on ',pick]);
text(50, 100, ['based on ',pick]);

%%
close (w)
toc
