% Sparse Group Lasso code
% The code for Connected Component Labeling is not provided here, so you can implement it yourself.
% ================================ Reference =============================
% Electrical Impedance Tomography for Thorax Image Using Connected Component Labeling and Sparse Group Lasso.
% In IEEE Transactions on Instrumentation and Measurement.
% ============= Author =============
% revised by Zekun Chen
%% base model 
clear all;
close all;
addpath('src');

load pha_data.mat
homg_data_norm=pha_data(:,1);%1-ref 
homg_data_norm=homg_data_norm/max(homg_data_norm(:));

inh_data_norm=pha_data(:,2);%2-data
inh_data_norm=inh_data_norm/max(inh_data_norm(:));

load logic.mat
homg_data_norm=homg_data_norm.*logic(:);
homg_data_norm(homg_data_norm==0)=[];
inh_data_norm=inh_data_norm.*logic(:);
inh_data_norm(inh_data_norm==0)=[];

load pha_data_sgt_n.mat
load .\src\S1.mat
load .\src\MM1.mat
load .\src\tranp1.mat

PixNum=64;%64*64
blocknum=3228;%block
range_s1=0.1;
range_s2=1;

dv = (-inh_data_norm+homg_data_norm)./max(homg_data_norm(:));
J = S;
%% NOSER+laplace
mask=zeros(blocknum,PixNum*PixNum);%mask
temp=tranp(:);
indx=find(temp==1);
for i=1:blocknum
    mask(i,indx(i))=1;
end

tranp(find(tranp==0))=nan;
hp  = 0.3;
d_filt = [0,1,0;1,-4,1;0,1,0];
Dl = convmtx2(d_filt, PixNum,PixNum) * mask';
Reg = Dl'*Dl;
RM  = left_divide((J'*J +  hp^2*Reg'*Reg),J');% left_divide
sol5=RM * dv;
tranp(tranp==0)=nan;
[sol5,PS] = mapminmax(sol5',range_s1,range_s2);
%% lasso model
% Parameter initialization
max_iter = 5;    
vareps    = 1e-4;
alpha    = 5;    
eta1=    0.001;  
eta2=    0.001;  
rho1    = 0.025;
rho2    = 0.25;
theta    = 0.2;
J = S;        
F=Reg;        
z=ones(blocknum,1); 
mu1=ones(blocknum,1); 
mu2=ones(208,1);   
mu1_pre=zeros(blocknum,1);
mu2_pre=zeros(208,1);
dv = -inh_data_norm+homg_data_norm; 
sigma_group=sgt(:,2:4);            % groundtruth
sgtt=sgt(:,1);
MM(:,3)=sgtt; 
ZZg=cal_vec2array(MM,blocknum,PixNum);
sigma=sol5';                
sigma_pre=zeros(blocknum,1); 
[~,N]=size(sigma_group);
w=ones(1,N);                       % group weight
%% ADMM solver
tic
sigma=ADMM_base(max_iter,vareps,sigma,mu1,mu2,alpha,eta1,eta2,rho1,rho2,theta,N,w,z,J,F,sigma_group,dv,blocknum);
toc
[sol,PS] = mapminmax(sigma',range_s1,range_s2);
MM(:,3)=sol;
ZZ1=cal_vec2array(MM,blocknum,PixNum);
ZZZ1=ZZ1.*tranp;
figure('color','white');
imagesc(ZZZ1,'AlphaData',tranp);
myColorMap=load('.\src\CM2.txt');
colormap(myColorMap);
axis square
axis off
% cb=colorbar;
% cb.Label.String = 'mS/cm';

S_rmse1=sqrt(immse(sol',sgtt))
S_2=ssim(sol',sgtt)
P_2=psnr(sol',sgtt) 
T2 = Tenengrad(ZZ1)
icc2=ICC(sol',sgtt)
M1=multissim(ZZ1,ZZg)

%% Reshape
load .\src\tranp256.mat
ZZd=ZZ1;
ZZd(ZZd(:)==0)=0.8;
ZZ2=imresize(ZZd,[256,256],'bilinear');

ZZ2=ZZ2.*tranp1;
figure('color','white');
imagesc(ZZ2,'AlphaData',tranp1);
myColorMap=load('.\src\CM2.txt');
colormap(myColorMap);
axis square
axis off
% cb=colorbar;
% cb.Label.String = 'mS/cm';