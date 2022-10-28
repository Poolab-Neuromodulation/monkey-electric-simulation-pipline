function simu_monkey(dn, subj, Ax, Ay, xl, xh, y, z, current)
% dn = '/home/poolab/code/monkey_simu_pipline/TEST/';
% subj = 'PP3';
tag = 'mesh';
cd(dn);
fn.label.tissue = fullfile(dn, 'label_tissue.nii');
current = str2num(current);
ind_e = reshape(find(current), [], length(find(current)));

segTouchup(fn.label.tissue);

coor_A = [str2num(Ax) str2num(Ay)];
boundary = [str2num(xl) str2num(xh) str2num(y) str2num(z)];
placeElec(fn.label.tissue,coor_A,boundary);

loc = {'Pz','Cz','Fz','Fpz',...
    'O1','C3','F1',...
    'O2','C4','F2',...
    'P1','P2',...
    'T1','T2'}';

labelParcel(subj);
[elem,node] = meshing(subj,tag);

prepareForGetDP(subj,node,elem,loc,tag);

load('conduct.mat')

tagloc = 'elec';

solveByGetDP(subj,current,conduct,1:length(loc),tag,tagloc,'');

visualizeResult(fn.label.tissue, subj, tag, tagloc, ind_e, loc);

waitfor(gcf);
end


function segTouchup(T1)

if nargin < 3 || isempty(isSmooth)
    isSmooth = 1; % smooth the tissue by default
end

if nargin < 4
    conn = [18 18 18 18 18 6];
end

disp('loading data...')
[dn,fn] = fileparts(T1);

nii.label.tissue = load_untouch_nii(fullfile(T1));

for i = 1:9
    mask(:,:,:,i) = nii.label.tissue.img == i;
end

[ex,mask] = binaryMaskGenerate(mask);
ex = removeFloatting(ex,18,'number',1);

% Fix csf continuity
mask = fixCSFcontunity(mask);

%% remove the floating objects
disp('removing disconnected voxels...')
thres = [1 1 3 1 500 8];
for i = 1:6
    if i == 4 || i ==6
        metho = 'size';
    else
        metho = 'number';
    end
    mask(:,:,:,i) = removeFloatting(mask(:,:,:,i),conn(i),metho,thres(i));
end
mask(:,:,:,6) = ex==1 | mask(:,:,:,6);


%% Generate unassigned voxels (empty voxels)
disp('generating and labeling empty voxels...')
empt = binaryMaskGenerate(mask);

itera = 0;
while any(empt(:))
    temp = uint8(mask)*255;
    sigma = 1;
    for i = 1:size(mask,4)
        temp(:,:,:,i) = imgaussfilt3(temp(:,:,:,i), sigma,'FilterSize',5);
    end
    prio = [9 8 7 6 2 1 3 4 5];
    [~,temp(:,:,:,prio)] = binaryMaskGenerate(temp(:,:,:,prio));
    mask = (repmat(empt,[1 1 1 9]) & temp) | mask;
    empt_novel = ~ (empt | any(mask,4));
    if sum(empt_novel,'all') == 0
        break;
    end
    empt = ~any(mask,4);
    itera = itera + 1;
    fprintf('   iteration %2d: %d empty volxels resting, %d novel empty volxels\n', ...
        itera,sum(empt(:)),sum(empt_novel(:)));
end
% Relabel each empty voxel to its nearest tissue type
% The Gaussian filter is used to calculate distances, and max operation
% relabels each empty voxel based on the distances.

disp('removing outside air...')
mask(:,:,:,6) = mask(:,:,:,6) & (~ex == 1);

% assign labels to tissues in this order: white,gray,csf,bone,skin,air
% note white-gray order differs from SPM outputs. Changed after ROAST V2.1
allMask = uint8(zeros(size(empt)));
for i = 1:size(mask,4)
    allMask(mask(:,:,:,i)) = i;
end

% save out
disp('Saving...');
nii.mask = nii.label.tissue;
nii.mask.img = uint8(allMask);
nii.mask.hdr.dime.scl_slope=1; % so that display of NIFTI will not alter the data
% In SPM results, this is 1/255, then all uint8 data will be displayed
% in the range of [0 1] % ANDY 2018-06-04
nii.mask.fileprefix = fullfile(dn,'masks');
nii.mask.hdr.hist.descrip = 'tissue masks';
save_untouch_nii(nii.mask,fullfile(dn,'masks.nii'));
end

function solveByGetDP(P,current,sigma,indUse,uniTag,uniTag2,LFtag)
% solveByGetDP(P,current,sigma,indUse,uniTag,LFtag)
% 
% Solve in getDP, a free FEM solver available at 
% http://getdp.info/
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017
% August 2019 adding lead field

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

load([dirname filesep baseFilename '_' uniTag '_usedElecArea.mat'],'area_elecNeeded');

numOfTissue = 6; % hard coded across ROAST.
numOfElec = length(area_elecNeeded);

fid = fopen([dirname filesep baseFilename '_' uniTag2 '.pro'],'w');

fprintf(fid,'%s\n\n','Group {');
fprintf(fid,'%s\n','white = Region[1];');
fprintf(fid,'%s\n','gray = Region[2];');
fprintf(fid,'%s\n','csf = Region[3];');
fprintf(fid,'%s\n','bone = Region[4];');
fprintf(fid,'%s\n','skin = Region[5];');
fprintf(fid,'%s\n','air = Region[6];');
% fprintf(fid,'%s\n','gel = Region[7];');
% fprintf(fid,'%s\n','elec = Region[8];');
for i=1:length(indUse)
    fprintf(fid,'%s\n',['gel' num2str(i) ' = Region[' num2str(numOfTissue+indUse(i)) '];']);
end
for i=1:length(indUse)
    fprintf(fid,'%s\n',['elec' num2str(i) ' = Region[' num2str(numOfTissue+numOfElec+indUse(i)) '];']);
end

gelStr = [];
elecStr = [];
usedElecStr = [];
for i=1:length(indUse)
%     fprintf(fid,'%s\n',['usedElec' num2str(i) ' = Region[' num2str(8+i) '];']);
    usedElecStr = [usedElecStr 'usedElec' num2str(i) ', '];
    fprintf(fid,'%s\n',['usedElec' num2str(i) ' = Region[' num2str(numOfTissue+2*numOfElec+indUse(i)) '];']);
    gelStr = [gelStr 'gel' num2str(i) ', '];
    elecStr = [elecStr 'elec' num2str(i) ', '];
end

% fprintf(fid,'%s\n','DomainC = Region[{white, gray, csf, bone, skin, air, gel, elec}];');
% fprintf(fid,'%s\n\n',['AllDomain = Region[{white, gray, csf, bone, skin, air, gel, elec, ' usedElecStr(1:end-2) '}];']);
fprintf(fid,'%s\n',['DomainC = Region[{white, gray, csf, bone, skin, air, ' gelStr elecStr(1:end-2) '}];']);
fprintf(fid,'%s\n\n',['AllDomain = Region[{white, gray, csf, bone, skin, air, ' gelStr elecStr usedElecStr(1:end-2) '}];']);
fprintf(fid,'%s\n\n','}');

fprintf(fid,'%s\n\n','Function {');
fprintf(fid,'%s\n',['sigma[white] = ' num2str(sigma.white) ';']);
fprintf(fid,'%s\n',['sigma[gray] = ' num2str(sigma.gray) ';']);
fprintf(fid,'%s\n',['sigma[csf] = ' num2str(sigma.csf) ';']);
fprintf(fid,'%s\n',['sigma[bone] = ' num2str(sigma.bone) ';']);
fprintf(fid,'%s\n',['sigma[skin] = ' num2str(sigma.skin) ';']);
fprintf(fid,'%s\n',['sigma[air] = ' num2str(sigma.air) ';']);
% fprintf(fid,'%s\n','sigma[gel] = 0.3;');
% fprintf(fid,'%s\n','sigma[elec] = 5.9e7;');
for i=1:length(indUse)
    fprintf(fid,'%s\n',['sigma[gel' num2str(i) '] = ' num2str(sigma.gel(indUse(i))) ';']);
end
for i=1:length(indUse)
    fprintf(fid,'%s\n',['sigma[elec' num2str(i) '] = ' num2str(sigma.electrode(indUse(i))) ';']);
end

for i=1:length(indUse)
    fprintf(fid,'%s\n',['du_dn' num2str(i) '[] = ' num2str(1000*current(indUse(i))/area_elecNeeded(indUse(i))) ';']);
end

fprintf(fid,'%s\n\n','}');

% fprintf(fid,'%s\n\n','Constraint {');
% fprintf(fid,'%s\n','{ Name ElectricScalarPotential; Type Assign;');
% fprintf(fid,'%s\n','  Case {');
% fprintf(fid,'%s\n','    { Region cathode; Value 0; }');
% fprintf(fid,'%s\n','  }');
% fprintf(fid,'%s\n\n','}');
% fprintf(fid,'%s\n\n','}');

fprintf(fid,'%s\n','Jacobian {');
fprintf(fid,'%s\n','  { Name Vol ;');
fprintf(fid,'%s\n','    Case {');
fprintf(fid,'%s\n','      { Region All ; Jacobian Vol ; }');
fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','  }');
fprintf(fid,'%s\n','  { Name Sur ;');
fprintf(fid,'%s\n','    Case {');
fprintf(fid,'%s\n','      { Region All ; Jacobian Sur ; }');
fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','  }');
fprintf(fid,'%s\n\n','}');

fprintf(fid,'%s\n','Integration {');
fprintf(fid,'%s\n','  { Name GradGrad ;');
fprintf(fid,'%s\n','    Case { {Type Gauss ;');
fprintf(fid,'%s\n','            Case { { GeoElement Triangle    ; NumberOfPoints  3 ; }');
fprintf(fid,'%s\n','                   { GeoElement Quadrangle  ; NumberOfPoints  4 ; }');
fprintf(fid,'%s\n','                   { GeoElement Tetrahedron ; NumberOfPoints  4 ; }');
fprintf(fid,'%s\n','                   { GeoElement Hexahedron  ; NumberOfPoints  6 ; }');
fprintf(fid,'%s\n','                   { GeoElement Prism       ; NumberOfPoints  9 ; } }');
fprintf(fid,'%s\n','           }');
fprintf(fid,'%s\n','         }');
fprintf(fid,'%s\n','  }');
fprintf(fid,'%s\n\n','}');

fprintf(fid,'%s\n','FunctionSpace {');
fprintf(fid,'%s\n','  { Name Hgrad_v_Ele; Type Form0;');
fprintf(fid,'%s\n','    BasisFunction {');
fprintf(fid,'%s\n','      // v = v  s   ,  for all nodes');
fprintf(fid,'%s\n','      //      n  n');
fprintf(fid,'%s\n','      { Name sn; NameOfCoef vn; Function BF_Node;');
fprintf(fid,'%s\n','        Support AllDomain; Entity NodesOf[ All ]; }');
fprintf(fid,'%s\n','    }');
% fprintf(fid,'%s\n','    Constraint {');
% fprintf(fid,'%s\n','      { NameOfCoef vn; EntityType NodesOf; ');
% fprintf(fid,'%s\n','        NameOfConstraint ElectricScalarPotential; }');
% fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','  }');
fprintf(fid,'%s\n\n','}');

fprintf(fid,'%s\n','Formulation {');
fprintf(fid,'%s\n','  { Name Electrostatics_v; Type FemEquation;');
fprintf(fid,'%s\n','    Quantity {');
fprintf(fid,'%s\n','      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }');
fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','    Equation {');
fprintf(fid,'%s\n','      Galerkin { [ sigma[] * Dof{d v} , {d v} ]; In DomainC; ');
fprintf(fid,'%s\n\n','                 Jacobian Vol; Integration GradGrad; }');

for i=1:length(indUse)
    
    fprintf(fid,'%s\n',['      Galerkin{ [ -du_dn' num2str(i) '[], {v} ]; In usedElec' num2str(i) ';']);
    fprintf(fid,'%s\n','                 Jacobian Sur; Integration GradGrad;}');
    
end

fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','  }');
fprintf(fid,'%s\n\n','}');

fprintf(fid,'%s\n','Resolution {');
fprintf(fid,'%s\n','  { Name EleSta_v;');
fprintf(fid,'%s\n','    System {');
fprintf(fid,'%s\n','      { Name Sys_Ele; NameOfFormulation Electrostatics_v; }');
fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','    Operation { ');
fprintf(fid,'%s\n','      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];');
fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','  }');
fprintf(fid,'%s\n\n','}');

fprintf(fid,'%s\n','PostProcessing {');
fprintf(fid,'%s\n','  { Name EleSta_v; NameOfFormulation Electrostatics_v;');
fprintf(fid,'%s\n','    Quantity {');
fprintf(fid,'%s\n','      { Name v; ');
fprintf(fid,'%s\n','        Value { ');
fprintf(fid,'%s\n','          Local { [ {v} ]; In AllDomain; Jacobian Vol; } ');
fprintf(fid,'%s\n','        }');
fprintf(fid,'%s\n','      }');
fprintf(fid,'%s\n','      { Name e; ');
fprintf(fid,'%s\n','        Value { ');
fprintf(fid,'%s\n','          Local { [ -{d v} ]; In AllDomain; Jacobian Vol; }');
fprintf(fid,'%s\n','        }');
fprintf(fid,'%s\n','      }');
fprintf(fid,'%s\n','    }');
fprintf(fid,'%s\n','  }');
fprintf(fid,'%s\n','}');

fprintf(fid,'%s\n\n','PostOperation {');
fprintf(fid,'%s\n','{ Name Map; NameOfPostProcessing EleSta_v;');
fprintf(fid,'%s\n','   Operation {');
if isempty(LFtag)
    fprintf(fid,'%s\n',['     Print [ v, OnElementsOf DomainC, File "' baseFilename '_' uniTag2 '_v.pos", Format NodeTable ];']);
end
fprintf(fid,'%s\n',['     Print [ e, OnElementsOf DomainC, Smoothing, File "' baseFilename '_' uniTag2 '_e' LFtag '.pos", Format NodeTable ];']);
fprintf(fid,'%s\n','   }');
fprintf(fid,'%s\n\n','}');
fprintf(fid,'%s\n','}');

fclose(fid);

str = computer('arch');
switch str
     case 'win64'
         solverPath = which('getdp.exe');
     case 'glnxa64'
         solverPath = which('getdp.exe');
         solverPath = solverPath(1:end-4);
     case 'maci64'
         solverPath = which('getdp.exe');
         solverPath = [solverPath(1:end-4), 'Mac'];
     otherwise
        error('Unsupported operating system!');
end

% cmd = [fileparts(which(mfilename)) filesep solverPath ' '...
%     fileparts(which(mfilename)) filesep dirname filesep baseFilename '_' uniTag '.pro -solve EleSta_v -msh '...
%     fileparts(which(mfilename)) filesep dirname filesep baseFilename '_' uniTag '_ready.msh -pos Map'];
cmd = [solverPath ' "' dirname filesep baseFilename '_' uniTag2 '.pro" -solve EleSta_v -msh "' dirname filesep baseFilename '_' uniTag '_ready.msh" -pos Map'];
try
    status = system(cmd);
catch
end

if status
    error('getDP solver cannot work properly on your system. Please check any error message you got.');
else % after solving, delete intermediate files
    delete([dirname filesep baseFilename '_' uniTag '.pre']);
    delete([dirname filesep baseFilename '_' uniTag '.res']);
end
end

function labelParcel(subj)
fn = 'masks_elec.nii'; 
fn_out = 'masks_final.nii';

nii = load_untouch_nii(fn);
C = load_untouch_nii(['CHARM_2_in_' subj '_padding.nii']);
S = load_untouch_nii(['SARM_2_in_' subj '_padding.nii']);
label_C = unique(C.img); label_C(1) = [];
label_S = unique(S.img); label_S(1) = [];

for i = 1:length(label_C)
    nii.img(nii.img == 2 & C.img == label_C(i)) = label_C(i) + 1000;
end

for i = 1:length(label_S)
    nii.img(nii.img == 2 & S.img == label_S(i)) = label_S(i) + 2000;
end
save_untouch_nii(nii,fn_out);

end

function [empt,mask] = binaryMaskGenerate(mask)
if islogical(mask)
    empt = ~any(mask,4);
    mask(:,:,:,1) = mask(:,:,:,1);
    for i = 2:size(mask,4)
        mask(:,:,:,i) = mask(:,:,:,i) & (~any(mask(:,:,:,1:i-1),4));
    end
else
    empt = mask(:,:,:,1);
    empt(:) = 0;
    data4D = cat(4,empt,mask);
    [~,ind] = max(data4D,[],4);
    empt = ind == 1;
    for i = 1:size(data4D,4)-1
         mask(:,:,:,i) = (ind == i+1);
    end
end 
end

function mask = removeFloatting(mask,conn,metho,thres)
arguments
    mask
    conn
    metho {mustBeMember(metho,{'size','number'})}
    thres
end
switch metho
    case 'size'
        mask = bwareaopen(mask,thres,conn);
    case 'number'
        siz = sizeOfObject(mask,conn);
        if length(siz) > thres
            thres = siz(thres+1)+1;
            mask = bwareaopen(mask,thres,conn);
        end
end
end

function mask = fixCSFcontunity(mask)
disp('fixing CSF continuity...')
se=ones(3,3,3);
empt = ~any(mask,4);
dcsf=imdilate(mask(:,:,:,3), se);
dbone=imdilate(mask(:,:,:,4), se);
contin=(empt&dcsf)|(dbone&mask(:,:,:,2));
mask(:,:,:,3)=mask(:,:,:,3)|contin;
mask_temp = mask(:,:,:,[3 4 1]);
[~,mask_temp] = binaryMaskGenerate(mask_temp);
mask(:,:,:,[3 4 1]) = mask_temp;
end

function loc = placeElec(mask,coor_A,boundary)

dn = fileparts(mask);
nii = load_untouch_nii(fullfile(dn,'masks.nii'));
M1 = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];
origin = round(inv(M1)*[0;0;0;1]+1);
i_o = origin(1);


scalp = any(nii.img,4);
scalp_surface = extrasurf(scalp,8);
ind = sub2ind(size(scalp),scalp_surface(:,1),scalp_surface(:,2),scalp_surface(:,3));
scalp_filled = false(size(scalp));
scalp_filled(ind) = true;
s_reverse = false(size(scalp));
s_reverse(boundary(1):boundary(2),boundary(3):end,boundary(4):end) = true;
scalp_surface_upper = scalp & s_reverse;

scalp_surface_upper = extrasurf(scalp_surface_upper,8);

l_c = scalp_surface_upper(scalp_surface_upper(:,1)==i_o,:);
l_c = sortedge(l_c);

num_elec_c = 5;

loc_c = perclocate(l_c,num_elec_c);

n = num_elec_c;
loc = loc_c;


while n > 1
    loc_stt = loc_c(1+(num_elec_c-n)/2,:);
    loc_end = loc_c(end-(num_elec_c-n)/2,:);
    v_ap = loc_stt - loc_end;
    v_norm = cross(v_ap,[1,0,0]); v_norm = v_norm/norm(v_norm);

    D = abs(v_norm*(scalp_surface_upper-loc_stt)')/norm(v_norm);
    D_thr = 0.55*0.5*3.^.5;

    % left elec
    ind_ll = find(abs(D)<=D_thr & scalp_surface_upper(:,1)'<=i_o);
    l_l = scalp_surface_upper(ind_ll,:); 
    l_l = sortedge(l_l,loc_stt); 
    loc_l = perclocate(l_l,n); 
    loc_l = loc_l(2:end-1,:);



    % right elec
    ind_lr = find(abs(D)<=D_thr & scalp_surface_upper(:,1)'>=i_o);
    l_r = scalp_surface_upper(ind_lr,:); 
    l_r = sortedge(l_r,loc_stt); 
    loc_r = perclocate(l_r,n); 
    loc_r = loc_r(2:end-1,:);


    loc = [loc;loc_l;loc_r];
    n = n-2;
end

ind_A = find(squeeze(scalp(:,coor_A(1),coor_A(2)))); ind_A = ind_A([1,end]);
loc_A = [ind_A [coor_A;coor_A]];
loc = [loc;loc_A];

loc(1,:) = []; % Occi lobe is removed

[~,indOnScalpSurf] = project2ClosestSurfacePoints(loc,scalp_surface,loc(floor(num_elec_c/2)+1,:));
elec_range = indOnScalpSurf(1:30,:)';
% [~,scalp_filled] = cleanScalp(any(nii.img,4),s_scalp);

% placing and model the electrodes
resolution = mean(nii.hdr.dime.pixdim(2:4));
% mean() here to handle anisotropic resolution; ugly. Maybe just
% resample MRI to isotropic in the very beginning?
[elec_C,gel_C] = placeAndModelElectrodes(loc,elec_range,scalp_surface,scalp_filled,resolution,origin(1:3));
mask_elec = cor2vol(elec_C,size(scalp));
mask_gel = cor2vol(gel_C,size(scalp));

label = unique(nii.img);
label(label<1 | label>=1000) = [];

for i = 1:length(loc)
    nii.img(nii.img==0 & mask_elec == i) = length(label) + length(loc) + i;
    nii.img(nii.img==0 & mask_gel == i) = length(label) + i;
end
save_untouch_nii(nii,'masks_elec.nii')


% for i = 1:size(loc,1)
%     scatter3(gel_C{i}(:,1),gel_C{i}(:,2),gel_C{i}(:,3),'filled','MarkerFaceColor',[1 0 0]); 
%     scatter3(elec_C{i}(:,1),elec_C{i}(:,2),elec_C{i}(:,3),'filled','MarkerFaceColor',[0 1 0]); 
% end
end

function l_s = sortedge(l,vargin)

if nargin>1
    p_stt = vargin;
    Dv = l - p_stt;
%     D = ;
    [~,ind] = min((Dv(:,1).^2 + Dv(:,2).^2 + Dv(:,3).^2).^5);
else
    p_stt = l(1,:); 
end

l_t = l;
l_s = nan(size(l));
l_s(1,:) = p_stt; 
l_t(1,:) = [];

for i = 2:length(l)
    [k,d] = dsearchn(l_t,l_s(i-1,:));
    l_s(i,:) = l_t(k,:);
    l_t(k,:) = [];
end

end

function [elec_allCoord,gel_allCoord] = placeAndModelElectrodes(elecLoc,elecRange,scalpCleanSurf,scalpFilled,res,origin)

disp('placing electrodes...')

% [Nx, Ny, Nz] = size(scalpFilled); % size of head in RAS orientation
scalpFilled(:,:,[1 end]) = 0; 
% scalpFilled(:,:,Nz) = 0; 
scalpFilled(:,[1 end],:) = 0; 
% scalpFilled(:,Ny,:) = 0; 
scalpFilled([1 end],:,:) = 0; 
% scalpFilled(Nx,:,:) = 0;


elec_allCoord = cell(size(elecLoc,1),1); gel_allCoord = cell(size(elecLoc,1),1);
% buffer for coordinates of each electrode and gel point
for i = 1:length(elecLoc) % size(elecLoc,1)    
    lcl = scalpCleanSurf(elecRange(i,:),:); % local scalp surface for each electrode

    [U,D] = eig(cov(lcl)); [~,ind] = min(diag(D));
    nv = U(:,ind)'; normal = nv/norm(nv); % Local normal for each electrode
    
    lenTry=1;
    testPointIn = round(elecLoc(i,:) - lenTry*normal);
    testPointOut = round(elecLoc(i,:) + lenTry*normal);
    while all(min([testPointIn;testPointOut])>0) && all(max([testPointIn;testPointOut])<=size(scalpFilled)) && ...
            ~xor(scalpFilled(testPointIn(1),testPointIn(2),testPointIn(3)),scalpFilled(testPointOut(1),testPointOut(2),testPointOut(3)))
        lenTry = lenTry+1;
        testPointIn = round(elecLoc(i,:) - lenTry*normal);
        testPointOut = round(elecLoc(i,:) + lenTry*normal);
    end
    if (elecLoc(i,:) - origin')*normal' < 0
        normal = -normal;
    end
    
    
            
    disc_radius = 6/res;
    disc_height = 2/res;
    
    dimTry = disc_radius;
    
    gel_out = elecLoc(i,:) +  2*disc_height*normal;
    electrode = gel_out + disc_height*normal;
    gel_in = gel_out - dimTry*normal; % coordinates of the boundaries of gel and electrode
    
    den = 2; % 2 points per pixel
    gel_coor = drawCylinder(0,disc_radius,gel_in,gel_out,den);
    elec_coor = drawCylinder(0,disc_radius,gel_out,electrode,den);
    % Use cylinders to model electrodes and gel, and calculate the coordinates of the points that make up the cylinder
    
    gel_coor = unique(round(gel_coor),'rows');
    elec_coor = unique(round(elec_coor),'rows'); % clean-up of the coordinates
    
    gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            
end
end

function loc = perclocate(l_c,n)

l_c = sortedge(l_c);

d = find(max(l_c,[],1)~=min(l_c,[],1));
if length(d)==2
    [bx,by,finalbreaks]=ncs2dapprox(l_c(:,d(1)),l_c(:,d(2)));
    % Approximation of 2-D Data by Natural Cubic Spline
    % http://www.mathworks.co.jp/matlabcentral/fileexchange/7617
    t = finalbreaks';
    l_c_b = nan(length(t),3);
    l_c_b(:,d) = [bx by];
    l_c_b(:,all(isnan(l_c_b),1)) = l_c(1,all(isnan(l_c_b),1));
else
    l_c_b = l_c;
    t = 1:length(l_c);
end
pp= spline(t,l_c_b');
range = linspace(1,t(end),10*t(end));
l_c_d = ppval(pp,range);
% distance_all = sum(sqrt(diff(l_c_d(1,:)).^2+diff(l_c_d(2,:)).^2));
dl_c_d = sqrt(diff(l_c_d(1,:)).^2+diff(l_c_d(2,:)).^2+diff(l_c_d(3,:)).^2);
d_l = [0 zeros(size(dl_c_d))];
for i = 2:length(d_l)
    d_l(i) = d_l(i-1) + dl_c_d(i-1);
end

frac = linspace(0,1,n);
[~,ind] = min(abs(d_l'/d_l(end)-frac),[],1);
loc = l_c_d(:,ind)';
end

function loc = removeloc(loc,loc_r,thr)

for i = 1:size(loc,1)
    for j = 1:size(loc_r,1)
        if norm(loc(i,:) - loc_r(j,:))<=thr
            loc(i,:) = NaN;
        end
    end
end
loc(all(isnan(loc),2),:) = [];
end

function vol = cor2vol(coor,siz)
vol = zeros(siz);
for i = 1:size(coor,1)
    vol(sub2ind(siz,coor{i}(:,1),coor{i}(:,2),coor{i}(:,3))) = i;
end
end

function s = extrasurf(v,vargin)
B = zeros([size(v) 3]);
    for i = 1:size(v,1)
        B(i,:,:,1) = bwperim(squeeze(v(i,:,:)),vargin(:));
    end
    for i = 1:size(v,2)
        B(:,i,:,2) = bwperim(squeeze(v(:,i,:)),vargin(:));
    end
    for i = 1:size(v,3)
        B(:,:,i,3) = bwperim(squeeze(v(:,:,i)),vargin(:));
    end    

    B = all(B>0,4);
    [s(:,1),s(:,2),s(:,3)] = ind2sub(size(v),find(B));
end

function [elem,node] = meshing(subj,tag)
% if isempty(dn), dn = pwd; end

% numOfTissue = 6; % hard coded across ROAST.  max(allMask(:));

fn = 'masks_final.nii';
nii = load_untouch_nii(fn);

allMask = uint8(zeros(size(nii.img)));

label = unique(nii.img); label(1) = [];

disp(length(label));
for i = 1:length(label)
    allMask(nii.img == label(i)) = i;
end

meshOpt = struct('radbound',3,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',5);

% opt.radbound = 5; % default 6, maximum surface element size
% opt.angbound = 30; % default 30, miminum angle of a surface triangle
% opt.distbound = 0.4; % default 0.5, maximum distance
% % between the center of the surface bounding circle and center of the element bounding sphere
% opt.reratio = 3; % default 3, maximum radius-edge ratio
% maxvol = 10; %100; % target maximum tetrahedral elem volume

[node,elem,face] = cgalv2m(allMask,meshOpt,meshOpt.maxvol);

% node1 = [node(:,1:3) ones(length(node),1)] * affine' + 0.5;
% node(:,1:3) = node1(:,1:3);
% node(:,1:3) = node(:,1:3) + 0.5; % then voxel space
% 
% for i=1:3, node(:,i) = node(:,i)*hdrInfo.pixdim(i); end
% Put mesh coordinates into pseudo-world space (voxel space but scaled properly
% using the scaling factors in the header) to avoid mistakes in
% solving. Putting coordinates into pure-world coordinates causes other
% complications. Units of coordinates are mm here. No need to convert into
% meter as voltage output from solver is mV.
% ANDY 2019-03-13

% for i = 1:3
%     node(:,i) = node(:,i)*data.hdr.dime.pixdim(1+i) + data.hdr.hist.(['qoffset_' 'x'+i-1]) + 0.5;
% end
% I have tired to transform the coordinates but failed. And I found I make
% a mistake for order of scaling and shift with 0.5 in voxel space
% lipzh@Shanghaitech.edu.cn on 2022-07-04

node(:,1:3) = node(:,1:3) + 0.5; % then voxel space

for i = 1:3
    node(:,i) = node(:,i)*nii.hdr.dime.pixdim(1+i);
end

save('mesh_parcell.mat','elem','face');

elem(elem(:,end)>length(label(label<1000)),end) = 2;
face(face(:,end)>length(label(label<1000)),end) = 2;

disp('saving mesh...')
% maskName = {'WHITE','GRAY','CSF','BONE','SKIN','AIR','STEEL','TI1','TI2','TI3'};

% maskName = cell(1,numOfTissue+numRecorder+1+numStimulator);
% maskName(1:numOfTissue) = {'WHITE','GRAY','CSF','BONE','SKIN','AIR'};
% for i=1:numRecorder, maskName{numOfTissue+i} = ['RECORDER' num2str(i)]; end
% maskName{numOfTissue+numRecorder+1} = 'CEMENT';
% for i=1:numStimulator, maskName{numOfTissue+numRecorder+1+i} = ['STIMULATOR' num2str(i)]; end
% savemsh(node(:,1:3),elem,[dn filesep 'mesh.msh'],maskName);
savemsh(node(:,1:3),elem,[subj '_' tag '.msh']);
save([subj '_' tag '.mat'],'node','elem','face');

% % visualize tissue by tissue
% figure;
% for i=1:length(unique(face(:,4)))
%     subplot(4,4,i)
%     title(maskName{i})
%     indElem = find(elem(:,5) == i);
%     indFace = find(face(:,4) == i);
%     plotmesh(node(:,1:3),face(indFace,:),elem(indElem,:),'LineStyle','none');
%     lightangle(-90,45)
%     lightangle(90,45)
%     lightangle(-90,-45)
%     rotate3d on;
%     view(-135,30);
% end
end

function prepareForGetDP(P,node,elem,elecNeeded,uniTag)
% prepareForGetDP(P,node,elem,elecNeeded,uniTag)
%
% Prepare to solve in getDP
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

% node = node + 0.5; already done right after mesh

% indNode_elecElm = elem(find(elem(:,5) == 8),1:4);
% X = zeros(size(indNode_elecElm,1),3);
% for e = 1:size(indNode_elecElm,1), X(e,:) = mean ( node(indNode_elecElm(e,:),1:3) ); end
% % figure; plot3(X(:,1),X(:,2),X(:,3),'r.');
% 
% Xt = round(X);
% label_elec = volume_elecLabel(sub2ind(size(volume_elecLabel),Xt(:,1),Xt(:,2),Xt(:,3)));
% indBad = find(label_elec==0);
% indGood = find(label_elec>0);
% [~,indOnGood] = map2Points(X(indBad,:),X(indGood,:),'closest'); % using nearest neighbor to fix bad labeling
% label_elec(indBad) = label_elec(indGood(indOnGood));
% 
% indNode_gelElm = elem(find(elem(:,5) == 7),1:4);
% X = zeros(size(indNode_gelElm,1),3);
% for e = 1:size(indNode_gelElm,1), X(e,:) = mean ( node(indNode_gelElm(e,:),1:3) ); end
% % figure; plot3(X(:,1),X(:,2),X(:,3),'r.');
% 
% Xt = round(X);
% label_gel = volume_gelLabel(sub2ind(size(volume_gelLabel),Xt(:,1),Xt(:,2),Xt(:,3)));
% indBad = find(label_gel==0);
% indGood = find(label_gel>0);
% [~,indOnGood] = map2Points(X(indBad,:),X(indGood,:),'closest'); % using nearest neighbor to fix bad labeling
% label_gel(indBad) = label_gel(indGood(indOnGood));
% 
% save([dirname filesep baseFilename '_' uniTag '_elecMeshLabels.mat'],'label_elec','label_gel');

numOfTissue = 6; % hard coded across ROAST.
numOfElec = length(elecNeeded);

element_elecNeeded = cell(numOfElec,1);
area_elecNeeded = zeros(numOfElec,1);

% resolution = mean(hdrInfo.pixdim);
% % mean() here to handle anisotropic rmesolution; ugly. Maybe just
% % resample MRI to isotropic in the very beginning?
% Not needed now when mesh coordinates are in pseudo-world space % ANDY 2019-03-13

warning('off','MATLAB:TriRep:PtsNotInTriWarnId');
for i=1:numOfElec
    
%     if isempty(indNode_elecElm(label_elec==i,:))
%         error(['Electrode ' elecNeeded{i} ' was not meshed properly. Reasons may be: 1) electrode size is too small so the mesher cannot capture it; 2) mesh resolution is not high enough. Consider using bigger electrodes or increasing the mesh resolution by specifying the mesh options.']);
%     end
%         
%     [faces_elec,verts_elec] = freeBoundary(TriRep(indNode_elecElm(label_elec==i,:),node(:,1:3)));
%     [faces_gel,verts_gel] = freeBoundary(TriRep(indNode_gelElm(label_gel==i,:),node(:,1:3)));
    
    indNode_gelElm = elem(find(elem(:,5) == numOfTissue+i),1:4);
    indNode_elecElm = elem(find(elem(:,5) == numOfTissue+numOfElec+i),1:4);
    
    if isempty(indNode_gelElm)
        error(['Gel under electrode ' elecNeeded{i} ' was not meshed properly. Reasons may be: 1) electrode size is too small so the mesher cannot capture it; 2) mesh resolution is not high enough. Consider using bigger electrodes or increasing the mesh resolution by specifying the mesh options.']);
    end
    
    if isempty(indNode_elecElm)
        error(['Electrode ' elecNeeded{i} ' was not meshed properly. Reasons may be: 1) electrode size is too small so the mesher cannot capture it; 2) mesh resolution is not high enough. Consider using bigger electrodes or increasing the mesh resolution by specifying the mesh options.']);
    end
    
    [faces_gel,verts_gel] = freeBoundary(TriRep(indNode_gelElm,node(:,1:3)));
    [faces_elec,verts_elec] = freeBoundary(TriRep(indNode_elecElm,node(:,1:3)));
    
    [~,iE,iG] = intersect(verts_elec,verts_gel,'rows');
    tempTag = ismember(faces_elec,iE);
    % faces_overlap = faces_elec(sum(tempTag,2)==3,:);
    faces_elecOuter = faces_elec(~(sum(tempTag,2)==3),:);
    [~,Loc] = ismember(verts_elec,node(:,1:3),'rows');
    element_elecNeeded{i} = Loc(faces_elecOuter);
    % calculate the surface area
    a = (verts_elec(faces_elecOuter(:, 2),:) - verts_elec(faces_elecOuter(:, 1),:)); %*resolution;
    b = (verts_elec(faces_elecOuter(:, 3),:) - verts_elec(faces_elecOuter(:, 1),:)); %*resolution;
    c = cross(a, b, 2);
    area_elecNeeded(i) = sum(0.5*sqrt(sum(c.^2, 2)));
    
end
if ~exist([dirname filesep baseFilename '_' uniTag '_usedElecArea.mat'],'file')
    save([dirname filesep baseFilename '_' uniTag '_usedElecArea.mat'],'area_elecNeeded');
end

if ~exist([dirname filesep baseFilename '_' uniTag '_ready.msh'],'file')
    
    disp('setting up boundary conditions...');
    
    fid_in = fopen([dirname filesep baseFilename '_' uniTag '.msh']);
    fid_out = fopen([dirname filesep baseFilename '_' uniTag '_ready.msh'],'w');
    
    numOfPart = length(unique(elem(:,5)));
    while ~feof(fid_in)
        s = fgetl(fid_in);
        
        if strcmp(s,'$Elements')
            fprintf(fid_out,'%s\n',s);
            s = fgetl(fid_in);
            numOfElem = str2num(s);
            fprintf(fid_out,'%s\n',num2str(numOfElem+size(cell2mat(element_elecNeeded),1)));
        elseif strcmp(s,'$EndElements')
            ii = 0;
            for j=1:numOfElec
                for i=1:size(element_elecNeeded{j},1)
                    
                    fprintf(fid_out,'%s \n',[num2str(numOfElem+i+ii) ' 2 2 ' num2str(numOfPart+j) ' ' num2str(numOfPart+j) ' ' num2str(element_elecNeeded{j}(i,1)) ' ' num2str(element_elecNeeded{j}(i,2)) ' ' num2str(element_elecNeeded{j}(i,3))]);
                    
                end
                ii = ii + i;
            end
            
            fprintf(fid_out,'%s\n',s);
        else
            fprintf(fid_out,'%s\n',s);
        end
    end
    
    fclose(fid_in);
    fclose(fid_out);
end   
end

function hm=plotmesh(node,varargin)
%
% hm=plotmesh(node,face,elem,opt)
%
% plot surface and volumetric meshes
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%      node: a node coordinate list, 3 columns for x/y/z; if node has a 
%            4th column, it will be used to set the color at each node.
%      face: a triangular surface face list; if face has a 4th column,
%            it will be used to separate the surface into 
%            sub-surfaces and display them in different colors;
%            face can be a cell array, each element of the array represents
%            a polyhedral facet of the mesh, if an element is an array with
%            two array subelements, the first one is the node index, the
%            second one is a scalar as the group id of the facet.
%      elem: a tetrahedral element list; if elem has a 5th column,
%            it will be used to separate the mesh into 
%            sub-domains and display them in different colors.
%      opt:  additional options for the plotting
%
%            for simple point plotting, opt can be markers
%            or color options, such as 'r.', or opt can be 
%            a logic statement to select a subset of the mesh,
%            such as 'x>0 & y+z<1', or an equation defining
%            a plane at which a mesh cross-section is plotted, for
%            example 'y=2*x'; opt can have more than one
%            items to combine these options, for example: 
%            plotmesh(...,'x>0','r.'); the range selector must
%            appear before the color/marker specifier
%
% in the event where all of the above inputs have extra settings related to 
% the color of the plot, the priorities are given in the following order:
%
%          opt > node(:,4) > elem(:,5) > face(:,4)
%
% output:
%   hm: handle or handles (vector) to the plotted surfaces
%
% example:
%
%   h=plotmesh(node,'r.');
%   h=plotmesh(node,'x<20','r.');
%   h=plotmesh(node,face);
%   h=plotmesh(node,face,'y>10');
%   h=plotmesh(node,face,'facecolor','r');
%   h=plotmesh(node,elem,'x<20');
%   h=plotmesh(node,elem,'x<20 & y>0');
%   h=plotmesh(node,face,elem);
%   h=plotmesh(node,face,elem,'linestyle','--');
%   h=plotmesh(node,elem,'z=20');
% 
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

selector=[];
opt=[];
face=[];
elem=[];

if(nargin>1)
   hasopt=0;
   for i=1:length(varargin)
   	if(ischar(varargin{i}))
		if(~isempty(regexp(varargin{i},'[x-zX-Z]')) && ~isempty(regexp(varargin{i},'[><=&|]')))
			selector=varargin{i};
			if(nargin>=i+1) opt=varargin(i+1:end); end
		else
			opt=varargin(i:end);
		end
		if(i==1)
			face=[];elem=[];
		elseif(i==2)
			if(iscell(varargin{1}) | size(varargin{1},2)<4)
				face=varargin{1}; elem=[];
			elseif(size(varargin{1},2)==4)
                faceid=unique(varargin{1}(:,4));
                if(length(faceid)==1)
                    face=varargin{1}; elem=[];
                elseif(any(hist(varargin{1}(:,4),unique(varargin{1}(:,4)))>50))
                    face=varargin{1}; elem=[];
                else
                    elem=varargin{1}; face=[];
                end
			else
				elem=varargin{1}; face=[];
			end
		elseif(i==3)
			face=varargin{1};
			elem=varargin{2};
		end
		hasopt=1;
		break;
	end
   end
   if(hasopt==0)
   	if(length(varargin)>=2)
		face=varargin{1};
		elem=varargin{2};
		if(length(varargin)>2) opt=varargin(3:end); end
	elseif(iscell(varargin{1}) | size(varargin{1},2)<4)
		face=varargin{1}; elem=[];
	elseif(size(varargin{1},2)==4)
	    faceid=unique(varargin{1}(:,4));
            if(length(faceid)==1)
	        face=varargin{1}; elem=[];
	    elseif(any(hist(varargin{1}(:,4),unique(varargin{1}(:,4)))>50))
                face=varargin{1}; elem=[];
	    else
                elem=varargin{1}; face=[];
	    end
	else
		elem=varargin{1}; face=[];
	end
   end
end

holdstate=ishold;
if(~holdstate)
    cla;
end
if(size(node,2)==4 && size(elem,2)==5)
    warning(['You have specified the node colors by both the 4th ' ...
            'and 5th columns of node and face inputs, respectively. ' ...
            'The node input takes priority']);
end
if(isempty(face) && isempty(elem))
   if(isempty(selector))
        if(isempty(opt))
   		h=plot3(node(:,1),node(:,2),node(:,3),'o');
	else
   		h=plot3(node(:,1),node(:,2),node(:,3),opt{:});
	end
   else
	x=node(:,1);
	y=node(:,2);
	z=node(:,3);
	idx=eval(['find(' selector ')']);
    if(~isempty(idx))
	    if(isempty(opt))
		h=plot3(node(idx,1),node(idx,2),node(idx,3),'o');
	    else
		h=plot3(node(idx,1),node(idx,2),node(idx,3),opt{:});
        end
    else
        warning('nothing to plot');
	end
   end
end

if(~isempty(face))
   hold on;
   if(isempty(selector))
        if(isempty(opt))
   		h=plotsurf(node,face);
	else
   		h=plotsurf(node,face,opt{:});
	end
   else
    if(iscell(face))
       cent=meshcentroid(node,face);
    else
       cent=meshcentroid(node,face(:,1:3));
    end
	x=cent(:,1);
    y=cent(:,2);
	z=cent(:,3);
    idx=eval(['find(' selector ')']);
    if(~isempty(idx))
        if(iscell(face))
            h=plotsurf(node,face(idx),opt{:});
        else
    		h=plotsurf(node,face(idx,:),opt{:});
        end
    else
        warning('no surface to plot');
	end
   end
end

if(~isempty(elem))
   hold on;
   if(isempty(selector))
        if(isempty(opt))
   		h=plottetra(node,elem);
	else
   		h=plottetra(node,elem,opt{:});
	end
   else
   cent=meshcentroid(node,elem(:,1:4));
   x=cent(:,1);
   y=cent(:,2);
   z=cent(:,3);
   if(regexp(selector,'='))
      if(size(node,2)==4)
          [cutpos,cutvalue,facedata]=qmeshcut(elem,node(:,1:3),node(:,4),selector);  
      elseif(size(node,2)==3)
          [cutpos,cutvalue,facedata]=qmeshcut(elem,node,node(:,3),selector);
      else
          error('plotmesh can only plot 3D tetrahedral meshes');
      end
      h=patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,'facecolor','interp',opt{:});
   else
      idx=eval(['find(' selector ')']);
      if(~isempty(idx))
	    if(isempty(opt))
		h=plottetra(node,elem(idx,:));
	    else
		h=plottetra(node,elem(idx,:),opt{:});
        end
      else
        warning('no tetrahedral element to plot');
	end
     end
   end
end

if(exist('h','var') & ~holdstate)
  hold off;
end
if(exist('h','var'))
  if(any(get(gca,'dataaspectratio')>1e8))
     view(3);
  end
  axis equal;
end
if(exist('h','var') & nargout>=1)
  hm=h;
end
end

function visualizeResult(P1, subj, tag, tagloc, ind_e, loc)
% visualizeRes(P1,P2,T2,node,elem,face,inCurrent,hdrInfo,uniTag,showAll,varargin)
%
% Display the simulation results. The 3D rendering is displayed in the
% world space, while the slice view is done in the voxel space.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018
% August 2019 callable by roast_target()

[dirname,baseFilename] = fileparts(P1);

data = load_untouch_nii([dirname filesep 'masks_final.nii']);
load(fullfile(dirname, [subj '_' tag '.mat']));
for i = 1:3
    node(:,i) = node(:,i) + data.hdr.hist.(['qoffset_' 'x'+i-1]);
end

indNode_grayFace = face(find(face(:,4) == 2),1:3);
indNode_grayElm = elem(find(elem(:,5) == 2),1:4);
    
fid = fopen([subj '_' tagloc '_e.pos']);
fgetl(fid);
C = textscan(fid,'%d %f %f %f');
fclose(fid);

% dataShow = [node(C{1},1:3), C_ef_mag];
color = nan(size(node,1),1);
color(C{1}) = sqrt(C{2}.^2+C{3}.^2+C{4}.^2);
dataShow = [node(:,1:3) color];

figure;
colormap(jet);
plotmesh(dataShow,indNode_grayFace,indNode_grayElm,'LineStyle','none');
hold on;
dataShowRange = [min(dataShow(unique(indNode_grayElm(:)),4)) prctile(dataShow(unique(indNode_grayElm(:)),4),95)];

axis off; rotate3d on;
caxis(dataShowRange);
lightangle(-90,45)
lightangle(90,45)
lightangle(-90,-45)

plotmesh(node(:,1:3),face(find(face(:,4) == 5),1:3),elem(find(elem(:,5) == 5),1:4),'LineStyle','none','FaceAlpha',0.01,'FaceColor',[.2 .2 .2])
for i = ind_e
    plotmesh(node(:,1:3),face(face(:,4)==i+length(loc)+6,1:3),elem(elem(:,5)==i+length(loc)+6,1:4),'LineStyle','none','FaceAlpha',.8,'FaceColor',[.8 .8 .8])
end

saveas(gcf, 'EF.jpg');

end