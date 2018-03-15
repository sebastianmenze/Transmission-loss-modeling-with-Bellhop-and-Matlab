%% make sim sources

addpath(genpath('C:\Users\Seb\Documents\passive_acoustic_work\matlab_functions'))
addpath(genpath('C:\Users\Seb\Documents\passive_acoustic_work\GridSphere'))
addpath(genpath('C:\Users\Seb\Documents\passive_acoustic_work\atWin10_2017_10\at_Compiled_W10_170915_64bit\Matlab\ReadWrite'))

%% read ETOPO


% latlim=[-72.5 -62];
% lonlim=[-50 15];
latlim=[-80 -45];
lonlim=[-65 25];

 etopo_location='C:\Users\Seb\Documents\passive_acoustic_work\ETOPO\etopo1_ice_c_f4.flt'
[Z_etopo, R_etopo] = etopo(etopo_location,1, latlim,lonlim);
 Z_etopo(Z_etopo>0)=0;
 
 %%
sim_sources.depth=10;

resolution_factor=6
resolution= 2 + ( 10 * (4 ^ resolution_factor) )

[sim_sources.lat,sim_sources.lon] = GridSphere(resolution) ;

% delete grid points ouside study area
ix_inside=sim_sources.lat>latlim(1) & sim_sources.lat<latlim(2) & sim_sources.lon>lonlim(1) & sim_sources.lon<lonlim(2) ;
sim_sources.lat(~ix_inside)=[];
sim_sources.lon(~ix_inside)=[];

% delete grid points on land
val = ltln2val(Z_etopo, R_etopo, sim_sources.lat,sim_sources.lon);
isOcean = val<0;

sim_sources.lat(~isOcean)=[];
sim_sources.lon(~isOcean)=[];

figure(1)
clf
set(gcf,'color','w')
hold on
worldmap(latlim,lonlim)
[coast]=load('coast');
plotm(coast.lat,coast.long,'-k')
plotm(sim_sources.lat,sim_sources.lon,'.k')

sim_sources.id=1:numel(sim_sources.lat);
sim_sources.id=sim_sources.id';
n_sources=numel(sim_sources.lat);
% get node spacing
for i=1:numel(sim_sources.lat)
    
 ix_without_i=1:numel(sim_sources.lat);
   
ix_without_i(ix_without_i==i)=[];
d_min(i)=min((distance(sim_sources.lat(i), sim_sources.lon(i),sim_sources.lat(ix_without_i), sim_sources.lon(ix_without_i))));
end

node_dist=mean(d_min);

%%  load global soundsspeed files

load('C:\Users\Seb\Documents\passive_acoustic_work\ssp_files\lev_latlonZ.mat')
load('C:\Users\Seb\Documents\passive_acoustic_work\ssp_files\lev_ann.mat')
c=c./100+1000;
lat=lat*0.1;
lon=lon*0.1;

for i_depth=1:numel(z)
[global_ssp(i_depth).z,global_ssp(i_depth).r]=geoloc2grid(lat, lon, c(:,i_depth), 1)
end
global_ssp_depth=z;
clear c z lat lon

%%
 S = referenceSphere('earth');
       angle_vector=0:5:180;

scenario_name='tl_bellhop'


%% loop through sim sources!
for i_source=1:n_sources
    
    
    write_folder=[scenario_name,'/raytracing_sim_sources_',num2str(i_source)];
    
if ~exist([scenario_name,'/slices_sim_source_',num2str(i_source),'.mat'])
       
for i_angle=1:numel(angle_vector)
    
s_lat=sim_sources.lat(i_source);
s_lon=sim_sources.lon(i_source);

[tlat1,tlon1] = track1(s_lat,s_lon,angle_vector(i_angle),50);
[tlat2,tlon2] = track1(s_lat,s_lon,rad2deg(deg2rad(angle_vector(i_angle))+pi),50);

% ixdel1= tlat1<latlim(2)+3 & tlat1>latlim(1)-3 & tlon1<lonlim(2)+3 & tlon1>lonlim(1)-3 ;
% ixdel2= tlat2<latlim(2)+3 & tlat2>latlim(1)-3 & tlon2<lonlim(2)+3 & tlon2>lonlim(1)-3 ;
% 
% tlat1(~ixdel1)=[];
% tlon1(~ixdel1)=[];
% tlat2(~ixdel2)=[];
% tlon2(~ixdel2)=[];

dis=distance(s_lat,s_lon,tlat1,tlon1);
[~,ixmd1]=max(dis);
dis=distance(s_lat,s_lon,tlat2,tlon2);
[~,ixmd2]=max(dis);

r_lat1=tlat1(ixmd1);
r_lon1=tlon1(ixmd1);
r_lat2=tlat2(ixmd2);
r_lon2=tlon2(ixmd2);

%%% for the ray boundaries:
[rtlat1,rtlon1] = track1(s_lat,s_lon,angle_vector(i_angle),60);
[rtlat2,rtlon2] = track1(s_lat,s_lon,rad2deg(deg2rad(angle_vector(i_angle))+pi),60);

ixdel1= rtlat1<latlim(2) & rtlat1>latlim(1) & rtlon1<lonlim(2) & rtlon1>lonlim(1) ;
ixdel2= rtlat2<latlim(2) & rtlat2>latlim(1) & rtlon2<lonlim(2) & rtlon2>lonlim(1) ;

rtlat1(~ixdel1)=[];
rtlon1(~ixdel1)=[];
rtlat2(~ixdel2)=[];
rtlon2(~ixdel2)=[];

dis=distance(s_lat,s_lon,rtlat1,rtlon1);
[~,ixmd1]=max(dis);
dis=distance(s_lat,s_lon,rtlat2,rtlon2);
[~,ixmd2]=max(dis);

rr_lat1=rtlat1(ixmd1);
rr_lon1=rtlon1(ixmd1);
rr_lat2=rtlat2(ixmd2);
rr_lon2=rtlon2(ixmd2);
%%%


% figure(1)
% clf
% set(gcf,'color','w')
% hold on
% worldmap('world')
% [coast]=load('coast');
% plotm(coast.lat,coast.long,'-k')
% plotm(sim_sources.lat,sim_sources.lon,'.k')
% plotm([s_lat],[s_lon],'or')
% 
% plotm(tlat1,tlon1,'.b')
% plotm(tlat2,tlon2,'.c')
% plotm([r_lat1],[r_lon1],'ob')
% plotm([r_lat2],[r_lon2],'oc')


%%%%% write bty file

clear zi ri profile_lat profile_lon
[zi,ri,~,~] = mapprofile(Z_etopo,R_etopo,[s_lat,r_lat1],[s_lon,r_lon1],'km');
is_nana=isnan(zi);
ri(is_nana)=[];
zi(is_nana)=[];
clear is_nana

% add negative range values
clear z2 r2
[z2,r2,~,~] = mapprofile(Z_etopo,R_etopo,[r_lat2,s_lat],[r_lon2,s_lon],'km');
is_nana=isnan(z2);
r2(is_nana)=[];
z2(is_nana)=[];
r2=r2-r2(end);
ri= [r2(1:end-1);ri] ;
zi= [z2(1:end-1);zi] ;


bty_name=['slice_',num2str(i_angle),'.bty'];

if exist(write_folder,'dir')
else
    mkdir(write_folder);
end

% Open file for writing
file_identifier = fopen([write_folder,'/',bty_name], 'w');

fprintf(file_identifier,'''L''\n');
fprintf(file_identifier,'%d\n',numel(ri));
for i=1:numel(ri)
    fprintf(file_identifier,'%6.1f %6.1f\n',ri(i),abs(zi(i)));
end
fclose(file_identifier);
clear file_identifier bty_name

deepest_point(i_angle)=abs(min(zi));

%% write ssp file

clear ssp_section
for i_depth=1:numel(global_ssp_depth)    
[ssp_section.ssp(i_depth,:),ssp_section.range,~,~] = mapprofile(global_ssp(i_depth).z,global_ssp(i_depth).r,[s_lat,r_lat1],[s_lon,r_lon1],'km');
end
% add negative range values
clear ssp2 r2
for i_depth=1:numel(global_ssp_depth)    
[ssp2(i_depth,:),r2,~,~] = mapprofile(global_ssp(i_depth).z,global_ssp(i_depth).r,[r_lat2,s_lat],[r_lon2,s_lon],'km');
end
r2=r2-r2(end);
ssp_section.ssp= [ssp2(:,1:end-1),ssp_section.ssp] ;
ssp_section.range=[r2(1:end-1); ssp_section.range];
    
% figure(2)
% clf
% subplot(211)
% plot(ri,zi,'-k')
% subplot(212)
% imagesc(ssp_section.range,global_ssp_depth,-ssp_section.ssp)
% colorbar

ssp_name=['slice_',num2str(i_angle),'.ssp'];

% Open file for writing
file_identifier = fopen([write_folder,'/',ssp_name], 'w');

fprintf(file_identifier,'%d\n',numel(ssp_section.range));
  row_string=[repmat('%6.2f ',1,numel(ssp_section.range)),'\n'];  
    fprintf(file_identifier,row_string,ssp_section.range);
for i_depth=1:numel(global_ssp_depth)
     fprintf(file_identifier,row_string,ssp_section.ssp(i_depth,:));  
   ssp_depth(i_angle,i_depth)=global_ssp_depth(i_depth);
   ssp_env(i_angle,i_depth)=ssp_section.ssp(i_depth,1)  ;
end
fclose(file_identifier);
clear file_identifier ssp_name



% replace woa end depth with deepest point
if ssp_depth(i_angle,end)<deepest_point(i_angle);
ssp_depth(i_angle,end)=deepest_point(i_angle);
end
%% run bellhop for 360deg slices
       
s_lat=sim_sources.lat(i_source);
s_lon=sim_sources.lon(i_source);
% r_lat=circle_lat(i_angle);
% r_lon=circle_lon(i_angle);

env_file=[write_folder,'/slice_',num2str(i_angle),'.env']
   disp(['run BELLHOP for: ' env_file])    
 
% Open file for writing
 frequency=20;
 sim_sources_depth=10;
 rec_depth=100;
 
file_identifier = fopen(env_file, 'w');

fprintf(file_identifier,'''%s''\n','Matlab generated env file'); % fclose(xfid)
fprintf(file_identifier,'%d ! FREQ(Hz)\n',frequency);
fprintf(file_identifier,'%d ! NMEDIA\n',1);
fprintf(file_identifier,'''%s'' \n','QVWT'); % for range dependent ssp
fprintf(file_identifier,'%d %6.1f %6.1f \n',[0,0,ssp_depth(i_angle,end)]);
for i_depth=1:numel(ssp_depth(i_angle,:))
fprintf(file_identifier,'%6.1f %6.1f /\n',[ssp_depth(i_angle,i_depth),ssp_env(i_angle,i_depth)]);
end
fprintf(file_identifier,'''A*'' 0.0\n');
fprintf(file_identifier,'%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n',[ssp_depth(i_angle,end),1800.0,0.0,2.0,0.1,0.0]);
fprintf(file_identifier,'%d \n',1);
fprintf(file_identifier,'%6.1f /\n',sim_sources_depth );
fprintf(file_identifier,'%d \n',1) ;
fprintf(file_identifier,'%6.1f /\n',rec_depth) ;
fprintf(file_identifier,'%d \n',1000) ;
fprintf(file_identifier,'%6.1f %6.1f /\n', [-deg2km(distance(s_lat,s_lon,rr_lat2,rr_lon2)),deg2km(distance(s_lat,s_lon,rr_lat1,rr_lon1))]) ;
fprintf(file_identifier,'''IB''\n');
% fprintf(file_identifier,'''IB''\n');

spacer=0.07;
angles=[-45:spacer:45,-135:spacer:-(180-spacer),135:spacer:180];

fprintf(file_identifier,'%d 0 \n',numel(angles)); % n rays choose by bellhop
% fprintf(file_identifier,'%d %d /\n',[-45 45]); % angle

fprintf(file_identifier,'%6.2f ',angles); % angle;
fprintf(file_identifier,'\n');

fprintf(file_identifier,'0.0 %6.1f %6.1f \n',[ssp_depth(i_angle,end), max(deg2km(distance(s_lat,s_lon,r_lat1,r_lon1)),deg2km(distance(s_lat,s_lon,r_lat2,r_lon2)))]);
fprintf(file_identifier,' \n');
fclose(file_identifier);

% RUN BELLHOP RAYTRACING MODEL
bellhop_location='C:\Users\Seb\Documents\passive_acoustic_work\atWin10_2017_10\at_Compiled_W10_170915_64bit\Bellhop\bellhop.exe';
system([ bellhop_location,' ',env_file(1:end-4)]);

 % get lat and lon of each of the 500 slice points
[slice_lat(i_angle,:),slice_lon(i_angle,:)] = track2(rr_lat2,rr_lon2,rr_lat1,rr_lon1,S,'degrees',1000);

[ PlotTitle, ~, freq, ~, Pos, pressure ] = read_shd( [env_file(1:end-4),'.shd'] );
slice_pressure(i_angle,:)=squeeze(pressure);
     
% figure(3)
% clf
% subplot(211)
% plot(ri,zi,'-k')
% subplot(212)
% plot(Pos.r.range,20*log10(squeeze(pressure)),'-k')

end
 
%%
figure(4)
clf
hold on
m_proj('lambert','long',lonlim,'lat',latlim);
m_gshhs_l('patch',[1 1 1]);
m_grid('xlabeldir','end','fontsize',10);

slice_tl=real(20*log10(slice_pressure));
slice_tl(slice_tl==-Inf)=NaN;
a=slice_tl(:);
b=slice_lat(:);
c=slice_lon(:);

plot_bins=100;
plot_value= a;
plot_range=[-150 -40];
[~,ind] = histc(plot_value,linspace(plot_range(1),plot_range(2),plot_bins));
cmap=parula(plot_bins);

for i=1:(plot_bins)    
  plot_cmap=cmap((i),:);   
m_plot(c(ind==i),b(ind==i),'.','markersize',20,'color',plot_cmap)
end
title(['Source ',num2str(i_source)])
%%
 print(gcf,'-dpng',[scenario_name,'\tl_map_source_',num2str(i_source)],'-r200')

 save([scenario_name,'/slices_sim_source_',num2str(i_source),'.mat'],'slice_lat','slice_lon','slice_pressure')
       clear slice_lat slice_lon slice_pressure slice_tl     

end 
end

  

