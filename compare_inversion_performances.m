
%% analyse test scenarios
addpath(genpath('C:\Users\Seb\Documents\passive_acoustic_work\'))
cd 'C:\Users\Seb\Documents\passive_acoustic_work\weddell_sea_scenarios'


latlim=[-80 -45];
lonlim=[-65 25];

latlim_nodes=[-72.5 -62];
lonlim_nodes=[-50 15];

sols=dir('solution_*.mat')

clear inversion 

[latgrid,longrid]=meshgrid(linspace(latlim_nodes(1),latlim_nodes(2),50),linspace(lonlim_nodes(1),lonlim_nodes(2),100));

for i_sol=1:numel(sols)

%     i_sol=find( strcmp( char(sols.name) ,'solution_000166.mat') )
    
solnume=str2num(sols(i_sol).name(end-9:end-4))
load(sols(i_sol).name);
load(['scenario',sols(i_sol).name(end-10:end)]);

ix_inside_nodes=true_sources.lat>latlim_nodes(1) & true_sources.lat<latlim_nodes(2) & true_sources.lon>lonlim_nodes(1) & true_sources.lon<lonlim_nodes(2);
inversion.inside_node_fraction(i_sol)=sum(true_sources.sl_p(ix_inside_nodes))/sum(true_sources.sl_p(:));

h=recorder.db_received;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));

inversion.rl_entropy(i_sol)=en_hist;
inversion.rl_entropy_value(i_sol)=en;

h=sim_sources.true_p;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.true_p_entropy(i_sol)=en_hist;
inversion.true_p_entropy_value(i_sol)=en;

%%%%%%%%%%% evaluation
[~,ix_best]=sort(likelihood);
ix_best=flip(ix_best);

% 2d correltion
F=scatteredInterpolant(sim_sources.lon,sim_sources.lat,sim_sources.true_p);
a1 = F(longrid,latgrid);
if sum(a1(:))==0
     a1(1,1)=1;
  %  a=rand(size(a));
end
F=scatteredInterpolant(true_sources.lon,true_sources.lat,true_sources.sl_p);
a2 = F(longrid,latgrid);
if sum(a2(:))==0
    a2(1,1)=1;
end
F=scatteredInterpolant(sim_sources.lon,sim_sources.lat,median(est_p(ix_best(1:3),:))');
b = F(longrid,latgrid);

inversion.r_smooth(i_sol) = corr2(a1,b);
inversion.r_detailed(i_sol) = corr2(a2,b);

% to find the trust zone
inversion.t(i_sol,:)=sim_sources.true_p;
inversion.s(i_sol,:)=median(est_p(ix_best(1:3),:))';


% figure(1)
% clf
% subplot(211)
%  imagesc(a2)
% subplot(212)
%  imagesc(b)
%  title(['smooth: ',num2str(inversion.r_smooth(i_sol)),' det: ',num2str(inversion.r_detailed(i_sol))])

% estimatied p entropy
h=est_p(ix_best(1),:);
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best(i_sol)=en_hist;
inversion.est_p_entropy_value_best(i_sol)=en;

h=median(est_p(ix_best(1:3),:));
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best3(i_sol)=en_hist;
inversion.est_p_entropy_value_best3(i_sol)=en;

h=hist(log(sse),10);
[x,y] = prepareCurveData(1:10,h);
[fitobject,gof] = fit(x,y,'exp1');
inversion.sse_fit(i_sol)=fitobject.b;

% tp=sim_sources.true_p;
% if sum(tp)==0
%     tp=tp+rand(size(tp))*.1
% end
% c=cov(sim_sources.true_p,est_p(ix_best(1),:));
% inversion.covar_best(i_sol)=c(1,2);
% 
% c=corrcoef(sim_sources.true_p,est_p(ix_best(1),:));
% inversion.p_r_best(i_sol)=c(1,2);

x=sim_sources.true_p;
if sum(x)==0
   y1=x; 
else
   y1=(x-min(x)) ./ (max(x)-min(x));   
end

x=est_p(ix_best(1),:);
y2=(x-min(x) )./ (max(x)-min(x));  
inversion.normsse_best(i_sol)=sum((y1'-y2).^2);

inversion.p_node_error(i_sol,:)=sim_sources.true_p' - est_p(ix_best(1),:);
inversion.p_node_error_norm(i_sol,:)=y1' - y2;

inversion.likelihood_best(i_sol)=likelihood(ix_best(1));
inversion.algorithm_sse_best(i_sol)=sse(ix_best(1));
%best 3
inversion.likelihood_best3(i_sol)=median(likelihood(ix_best(1:3)));
inversion.algorithm_sse_best3(i_sol)=median(sse(ix_best(1:3)));

% c=cov(sim_sources.true_p,median(est_p(ix_best(1:3),:)));
% inversion.covar_best3(i_sol)=c(1,2);
% 
% c=corrcoef(sim_sources.true_p,median(est_p(ix_best(1:3),:)));
% inversion.p_r_best3(i_sol)=c(1,2);

x=median(est_p(ix_best(1:3),:));
y2= (x-min(x)) ./ (max(x)-min(x));  
inversion.normsse_best3(i_sol)=sum((y1'-y2).^2);

% figure(1)
% clf
% subplot(121)
% box on
% hold on
% grid on
% plot(y1,y2,'.k')
% 
% subplot(122)
% box on
% hold on
% grid on
% plot(sim_sources.true_p,median(est_p(ix_best(1:3),:)),'.k')
% 
% pause

inversion.p_node_error_best3(i_sol,:)=sim_sources.true_p' - median(est_p(ix_best(1:3),:));
inversion.p_node_error_norm_best3(i_sol,:)=y1' - y2;

% indiv solutions
inversion.p_sse_best(i_sol)=sum( ( sim_sources.true_p' - est_p(ix_best(1),:) ).^2 );
for k=1:numel(likelihood)
    inversion.p_sse_solutions(i_sol,k)=sum( ( sim_sources.true_p' - est_p(k,:) ).^2 );
end

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% map plot

% figure(1)
% clf
% 
% subplot(221)
% hold on
% 
% m_proj('lambert','long',lonlim,'lat',latlim);
% 
% m_pcolor(ran.lon,ran.lat,ran.p);
% shading flat
% m_gshhs_l('patch',[1 1 1]);
% m_grid('xlabeldir','end','fontsize',10);
% 
%  plot_value=sim_sources.true_p;
% % plot_value=20*log10(sim_sources.true_p);
% % plot_value(plot_value==-Inf)=0;
% plot_bins=30;
%  plot_range=[min(plot_value),max(plot_value)];
% %  plot_range=[0,p_sum_smooth_max];
% jet_cmap=cmocean('matter',plot_bins);
% [~,ind] = histc(plot_value,linspace(plot_range(1),plot_range(2),plot_bins));
%  ind(ind==0)=plot_bins;
% 
% for j=1:plot_bins  
% plot_cmap=jet_cmap(j,:); 
% ix_color=ind==j;
% m_plot(sim_sources.lon(ix_color),sim_sources.lat(ix_color),'.','markersize',20,'color',plot_cmap)
% end
% 
% %%%% plot receivers
% plot_bins=100;
% plot_value= recorder.db_received;
% plot_range=[min(plot_value),max(plot_value)];
% [~,ind] = histc(plot_value,linspace(plot_range(1),plot_range(2),plot_bins));
% ind(ind==0)=plot_bins;
% cmap=cool(plot_bins);
% 
% for i=1:numel(plot_value)    
%   plot_cmap=cmap(ind(i),:);   
% m_plot(recorder.lon(i),recorder.lat(i),'.','markersize',20,'color',plot_cmap)
% end
% 
% 
% title(['True sum(p): ',num2str(sum(sim_sources.true_p),'%i')])
% 
% %%% plot best
% 
% subplot(222)
% hold on
% 
% m_proj('lambert','long',lonlim,'lat',latlim);
% 
% m_pcolor(ran.lon,ran.lat,ran.p);
% shading flat
% m_gshhs_l('patch',[1 1 1]);
% m_grid('xlabeldir','end','fontsize',10);
% 
% plot_value=est_p(ix_best(1),:);
% plot_bins=30;
%  plot_range=[min(plot_value),max(plot_value)];
% %  plot_range=[0,p_sum_smooth_max];
% jet_cmap=cmocean('matter',plot_bins);
% [~,ind] = histc(plot_value,linspace(plot_range(1),plot_range(2),plot_bins));
%  ind(ind==0)=plot_bins;
% 
% for j=1:plot_bins
% plot_cmap=jet_cmap(j,:); 
% ix_color=ind==j;
% m_plot(sim_sources.lon(ix_color),sim_sources.lat(ix_color),'.','markersize',20,'color',plot_cmap)
% end
% title(['Best Sol Norm SSE: ',num2str(inversion.normsse_best(i_sol))])
% 
% %%%%%%%%median sol
% 
% subplot(223)
% hold on
% 
% m_proj('lambert','long',lonlim,'lat',latlim);
% 
% m_pcolor(ran.lon,ran.lat,ran.p);
% shading flat
% m_gshhs_l('patch',[1 1 1]);
% m_grid('xlabeldir','end','fontsize',10);
% 
% plot_value=median(est_p(ix_best(1:3),:));
% plot_bins=30;
%  plot_range=[min(plot_value),max(plot_value)];
% %  plot_range=[0,p_sum_smooth_max];
% jet_cmap=cmocean('matter',plot_bins);
% [~,ind] = histc(plot_value,linspace(plot_range(1),plot_range(2),plot_bins));
%  ind(ind==0)=plot_bins;
% 
% for j=1:plot_bins
% plot_cmap=jet_cmap(j,:); 
% ix_color=ind==j;
% m_plot(sim_sources.lon(ix_color),sim_sources.lat(ix_color),'.','markersize',20,'color',plot_cmap)
% end
% title(['Median of best 3 Sol Norm SSE: ',num2str(inversion.normsse_best3(i_sol))])
% 
% subplot(224)
% hold on
% grid on
% bar(est_p_sum,sse,.4)
% xlabel('sum(p_{est})')
% ylabel('log(SSE)')
% set(gca,'yscale','log')
% title({['RL entropy: ',num2str(inversion.rl_entropy_value(i_sol))],['Node pressure distribution entropy: ',num2str(inversion.true_p_entropy(i_sol))],['Node pressure entropy: ',num2str(inversion.true_p_entropy_value(i_sol))],['Inside node fraction: ',num2str(inversion.inside_node_fraction(i_sol))]})
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% print(gcf,'-dpng',['inversion_map2_',sols(i_sol).name(end-10:end-4)],'-r200')
% %%%%%%%%%%%%%%%%%%%%%%%%


end


clear t s
for i=1:size(inversion.t,1)
    x=inversion.t(i,:);
    if sum(x)==0
     t(i,:)=x;   
    else
     x=(x-min(x) )./ (max(x)-min(x));  
    t(i,:)=x;
    end
    x=inversion.s(i,:);
    x=(x-min(x) )./ (max(x)-min(x));  
    s(i,:)=x;
end
for i=1:size(t,2)
c=corrcoef(t(:,i),s(:,i));
inversion.r_map_norm(i)=c(1,2);
end

for i=1:size(inversion.t,2)
c=corrcoef(inversion.t(:,i),inversion.s(:,i));
inversion.r_map(i)=c(1,2);
end

%%%%%%%% make map of trust zones

figure(4)
clf
hold on

m_proj('lambert','long',lonlim,'lat',latlim);
m_gshhs_l('patch',[1 1 1]);
m_grid('xlabeldir','end','fontsize',10);

plot_value=inversion.r_map;
plot_value(isnan(plot_value))=0;
plot_value(plot_value<0)=0;
plot_bins=30;
 plot_range=[0,1];
%  plot_range=[0,p_sum_smooth_max];
jet_cmap=cmocean('amplitude',plot_bins);
[~,ind] = histc(plot_value,linspace(plot_range(1),plot_range(2),plot_bins));
 ind(ind==0)=plot_bins;

for j=1:plot_bins  
plot_cmap=jet_cmap(j,:); 
ix_color=ind==j;
m_plot(sim_sources.lon(ix_color),sim_sources.lat(ix_color),'.','markersize',20,'color',plot_cmap)
end

inv{1}=inversion;



%% analyse test scenarios

clear inversion 

addpath(genpath('C:\Users\Seb\Documents\passive_acoustic_work\'))
cd 'C:\Users\Seb\Documents\passive_acoustic_work\weddell_sea_scenarios_2'


latlim=[-80 -45];
lonlim=[-65 25];

latlim_nodes=[-72.5 -62];
lonlim_nodes=[-50 15];
[latgrid,longrid]=meshgrid(linspace(latlim_nodes(1),latlim_nodes(2),50),linspace(lonlim_nodes(1),lonlim_nodes(2),100));

sols=dir('solution_*.mat')

for i_sol=1:numel(sols)

%     i_sol=find( strcmp( char(sols.name) ,'solution_000166.mat') )
    
solnume=str2num(sols(i_sol).name(end-9:end-4))
load(sols(i_sol).name);
load(['scenario',sols(i_sol).name(end-10:end)]);

ix_inside_nodes=true_sources.lat>latlim_nodes(1) & true_sources.lat<latlim_nodes(2) & true_sources.lon>lonlim_nodes(1) & true_sources.lon<lonlim_nodes(2);
inversion.inside_node_fraction(i_sol)=sum(true_sources.sl_p(ix_inside_nodes))/sum(true_sources.sl_p(:));

h=recorder.db_received;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));

inversion.rl_entropy(i_sol)=en_hist;
inversion.rl_entropy_value(i_sol)=en;

h=sim_sources.true_p;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.true_p_entropy(i_sol)=en_hist;
inversion.true_p_entropy_value(i_sol)=en;

%%%%%%%%%%% evaluation
[~,ix_best]=sort(likelihood);
ix_best=flip(ix_best);

% 2d correltion
F=scatteredInterpolant(sim_sources.lon,sim_sources.lat,sim_sources.true_p);
a1 = F(longrid,latgrid);
if sum(a1(:))==0
     a1(1,1)=1;
  %  a=rand(size(a));
end
F=scatteredInterpolant(true_sources.lon,true_sources.lat,true_sources.sl_p);
a2 = F(longrid,latgrid);
if sum(a2(:))==0
    a2(1,1)=1;
end
F=scatteredInterpolant(sim_sources.lon,sim_sources.lat,median(est_p(ix_best(1:3),:))');
b = F(longrid,latgrid);

inversion.r_smooth(i_sol) = corr2(a1,b);
inversion.r_detailed(i_sol) = corr2(a2,b);

% estimatied p entropy
h=est_p(ix_best(1),:);
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best(i_sol)=en_hist;
inversion.est_p_entropy_value_best(i_sol)=en;

h=median(est_p(ix_best(1:3),:));
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best3(i_sol)=en_hist;
inversion.est_p_entropy_value_best3(i_sol)=en;

h=hist(log(sse),10);
[x,y] = prepareCurveData(1:10,h);
[fitobject,gof] = fit(x,y,'exp1');
inversion.sse_fit(i_sol)=fitobject.b;

% tp=sim_sources.true_p;
% if sum(tp)==0
%     tp=tp+rand(size(tp))*.1
% end
% c=cov(sim_sources.true_p,est_p(ix_best(1),:));
% inversion.covar_best(i_sol)=c(1,2);
% 
% c=corrcoef(sim_sources.true_p,est_p(ix_best(1),:));
% inversion.p_r_best(i_sol)=c(1,2);

x=sim_sources.true_p;
if sum(x)==0
   y1=x; 
else
   y1=(x-min(x)) ./ (max(x)-min(x));   
end

x=est_p(ix_best(1),:);
y2=(x-min(x) )./ (max(x)-min(x));  
inversion.normsse_best(i_sol)=sum((y1'-y2).^2);

inversion.p_node_error(i_sol,:)=sim_sources.true_p' - est_p(ix_best(1),:);
inversion.p_node_error_norm(i_sol,:)=y1' - y2;

inversion.likelihood_best(i_sol)=likelihood(ix_best(1));
inversion.algorithm_sse_best(i_sol)=sse(ix_best(1));
%best 3
inversion.likelihood_best3(i_sol)=median(likelihood(ix_best(1:3)));
inversion.algorithm_sse_best3(i_sol)=median(sse(ix_best(1:3)));

% c=cov(sim_sources.true_p,median(est_p(ix_best(1:3),:)));
% inversion.covar_best3(i_sol)=c(1,2);
% 
% c=corrcoef(sim_sources.true_p,median(est_p(ix_best(1:3),:)));
% inversion.p_r_best3(i_sol)=c(1,2);

x=median(est_p(ix_best(1:3),:));
y2= (x-min(x)) ./ (max(x)-min(x));  
inversion.normsse_best3(i_sol)=sum((y1'-y2).^2);

inversion.p_node_error_best3(i_sol,:)=sim_sources.true_p' - median(est_p(ix_best(1:3),:));
inversion.p_node_error_norm_best3(i_sol,:)=y1' - y2;

% indiv solutions
inversion.p_sse_best(i_sol)=sum( ( sim_sources.true_p' - est_p(ix_best(1),:) ).^2 );
for k=1:numel(likelihood)
    inversion.p_sse_solutions(i_sol,k)=sum( ( sim_sources.true_p' - est_p(k,:) ).^2 );
end

end

inv{2}=inversion;

%% tests 3

clear inversion 

addpath(genpath('C:\Users\Seb\Documents\passive_acoustic_work\'))
cd 'C:\Users\Seb\Documents\passive_acoustic_work\weddell_sea_scenarios_3'


latlim=[-80 -45];
lonlim=[-65 25];

latlim_nodes=[-72.5 -62];
lonlim_nodes=[-50 15];
[latgrid,longrid]=meshgrid(linspace(latlim_nodes(1),latlim_nodes(2),50),linspace(lonlim_nodes(1),lonlim_nodes(2),100));

sols=dir('solution_*.mat')

for i_sol=1:numel(sols)

%     i_sol=find( strcmp( char(sols.name) ,'solution_000166.mat') )
    
solnume=str2num(sols(i_sol).name(end-9:end-4))
load(sols(i_sol).name);
load(['scenario',sols(i_sol).name(end-10:end)]);

ix_inside_nodes=true_sources.lat>latlim_nodes(1) & true_sources.lat<latlim_nodes(2) & true_sources.lon>lonlim_nodes(1) & true_sources.lon<lonlim_nodes(2);
inversion.inside_node_fraction(i_sol)=sum(true_sources.sl_p(ix_inside_nodes))/sum(true_sources.sl_p(:));

h=recorder.db_received;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));

inversion.rl_entropy(i_sol)=en_hist;
inversion.rl_entropy_value(i_sol)=en;

h=sim_sources.true_p;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.true_p_entropy(i_sol)=en_hist;
inversion.true_p_entropy_value(i_sol)=en;

%%%%%%%%%%% evaluation
[~,ix_best]=sort(likelihood);
ix_best=flip(ix_best);

%%%%%%%%%%%%%%%%%%%%% just inside grid
ix_inside_simnodes=sim_sources.lat>latlim_nodes(1) & sim_sources.lat<latlim_nodes(2) & sim_sources.lon>lonlim_nodes(1) & sim_sources.lon<lonlim_nodes(2);

% 2d correltion
F=scatteredInterpolant(sim_sources.lon,sim_sources.lat,sim_sources.true_p);
a1 = F(longrid,latgrid);
if sum(a1(:))==0
     a1(1,1)=1;
  %  a=rand(size(a));
end
F=scatteredInterpolant(true_sources.lon,true_sources.lat,true_sources.sl_p);
a2 = F(longrid,latgrid);
if sum(a2(:))==0
    a2(1,1)=1;
end
F=scatteredInterpolant(sim_sources.lon,sim_sources.lat,median(est_p(ix_best(1:3),:))');
b = F(longrid,latgrid);
if sum(b(:))==0
    b(1,1)=1;
end
inversion.r_smooth(i_sol) = corr2(a1,b);
inversion.r_detailed(i_sol) = corr2(a2,b);

% to find the trust zone
inversion.t(i_sol,:)=sim_sources.true_p;
inversion.s(i_sol,:)=median(est_p(ix_best(1:3),:))';

% estimatied p entropy
h=est_p(ix_best(1),ix_inside_simnodes);
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best(i_sol)=en_hist;
inversion.est_p_entropy_value_best(i_sol)=en;

h=median(est_p(ix_best(1:3),ix_inside_simnodes));
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best3(i_sol)=en_hist;
inversion.est_p_entropy_value_best3(i_sol)=en;

x=sim_sources.true_p(ix_inside_simnodes);
if sum(x)==0
   y1=x; 
else
   y1=(x-min(x)) ./ (max(x)-min(x));   
end

x=est_p(ix_best(1),ix_inside_simnodes);
y2=(x-min(x) )./ (max(x)-min(x));  
inversion.normsse_best(i_sol)=sum((y1'-y2).^2);

inversion.p_node_error(i_sol,:)=sim_sources.true_p(ix_inside_simnodes)' - est_p(ix_best(1),ix_inside_simnodes);
inversion.p_node_error_norm(i_sol,:)=y1' - y2;

x=median(est_p(ix_best(1:3),ix_inside_simnodes));
y2= (x-min(x)) ./ (max(x)-min(x));  
inversion.normsse_best3(i_sol)=sum((y1'-y2).^2);

% figure(1)
% clf
% subplot(121)
% box on
% hold on
% grid on
% plot(y1,y2,'.k')
% xlim([0,1])
% ylim([0,1])
% plot([0,1],[0,1],'-k')
% 
% subplot(122)
% box on
% hold on
% grid on
% plot(sim_sources.true_p(ix_inside_simnodes),median(est_p(ix_best(1:3),ix_inside_simnodes)),'.k')
% 
% pause

inversion.p_node_error_best3(i_sol,:)=sim_sources.true_p(ix_inside_simnodes)' - median(est_p(ix_best(1:3),ix_inside_simnodes));
inversion.p_node_error_norm_best3(i_sol,:)=y1' - y2;

% indiv solutions
inversion.p_sse_best(i_sol)=sum( ( sim_sources.true_p(ix_inside_simnodes)' - est_p(ix_best(1),ix_inside_simnodes) ).^2 );
for k=1:numel(likelihood)
    inversion.p_sse_solutions(i_sol,k)=sum( ( sim_sources.true_p(ix_inside_simnodes)' - est_p(k,ix_inside_simnodes) ).^2 );
end



end

% clear t s
% for i=1:size(inversion.t,1)
%     x=inversion.t(i,:);
%     if sum(x)==0
%      t(i,:)=x;   
%     else
%      x=(x-min(x) )./ (max(x)-min(x));  
%     t(i,:)=x;
%     end
%     x=inversion.s(i,:);
%     x=(x-min(x) )./ (max(x)-min(x));  
%     s(i,:)=x;
% end
% for i=1:size(t,2)
% c=corrcoef(t(:,i),s(:,i));
% inversion.r_map_norm(i)=c(1,2);
% end


%%%%%%%%%%%%%%%
inv{3}=inversion;

%% trust correlation map


for i=1:size(inv{3}.t,2)
c=corrcoef(inv{3}.t(:,i),inv{3}.s(:,i));
inv{3}.r_map(i)=c(1,2);
end

%%%%%%%% make map of trust zones

figure(4)
clf
hold on

m_proj('lambert','long',lonlim,'lat',latlim);
m_gshhs_l('patch',[1 1 1]);
m_grid('xlabeldir','end','fontsize',10);

plot_value=inv{3}.r_map;
plot_value(isnan(plot_value))=0;
plot_value(plot_value<0)=0;
plot_bins=30;
 plot_range=[0,1];
%  plot_range=[0,p_sum_smooth_max];
jet_cmap=cmocean('amplitude',plot_bins);
[~,ind] = histc(plot_value,linspace(plot_range(1),plot_range(2),plot_bins));
 ind(ind==0)=plot_bins;

for j=1:plot_bins  
plot_cmap=jet_cmap(j,:); 
ix_color=ind==j;
m_plot(sim_sources.lon(ix_color),sim_sources.lat(ix_color),'.','markersize',20,'color',plot_cmap)
end

m_plot(recorder.lon,recorder.lat,'.','markersize',20,'color','b')

% 
% imagefolder='C:\Users\Seb\Documents\passive_acoustic_work\figures\'
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',[imagefolder,'Map of scenario correlation factors (Algorithm trust area)'],'-r300') 

%% histogram norm sse

figure(2)
clf
hold on
grid on

bins=0:40;
b=histc(inv{1}.normsse_best3,bins);
plot(bins,b./sum(b),'o-r')
b=histc(inv{2}.normsse_best3,bins);
plot(bins,b./sum(b),'o-b')
b=histc(inv{3}.normsse_best3,bins);
plot(bins,b./sum(b),'o-c')

xlabel('Normalised SSE bewten true and estimated p')
ylabel('n')

%% histogram r

figure(2)
clf

subplot(211)
hold on
grid on

bins=0:.05:1;

b=histc(inv{1}.r_smooth,bins);
plot(bins,b./sum(b),'o-r')
b=histc(inv{2}.r_smooth,bins);
plot(bins,b./sum(b),'o-b')
b=histc(inv{3}.r_smooth,bins);
plot(bins,b./sum(b),'o-c')

xlabel('r')
ylabel('n')

subplot(212)
hold on
grid on

b=histc(inv{1}.r_detailed,bins);
plot(bins,b./sum(b),'o-r')
b=histc(inv{2}.r_detailed,bins);
plot(bins,b./sum(b),'o-b')
b=histc(inv{3}.r_detailed,bins);
plot(bins,b./sum(b),'o-c')

xlabel('r')
ylabel('n')

%%

figure(11)
clf

subplot(121)
hold on
box on
grid on

scatter(inv{1}.true_p_entropy,inv{1}.normsse_best,10,inv{1}.inside_node_fraction,'filled')
xlabel('Entropy of true node pressure distribution')
ylabel('Normalized SSE')
cb=colorbar
ylabel(cb,'Inside node fraction')


subplot(122)
hold on
box on
grid on

scatter(inv{2}.true_p_entropy,inv{2}.normsse_best,10,inv{2}.inside_node_fraction,'filled')
xlabel('Entropy of true node pressure distribution')
ylabel('Normalized SSE')
cb=colorbar
ylabel(cb,'Inside node fraction')

%%

figure(11)
clf

subplot(121)
hold on
box on
grid on

scatter(inv{1}.est_p_entropy_value_best3,inv{1}.normsse_best3,10,inv{1}.inside_node_fraction,'filled')

xlabel('Entropy of estimated node pressure')
ylabel('Normalised SSE')
cb=colorbar
ylabel(cb,'Inside node fraction')


subplot(122)
hold on
box on
grid on

scatter(inv{2}.est_p_entropy_value_best3,inv{2}.normsse_best3,10,inv{2}.inside_node_fraction,'filled')

xlabel('Entropy of estimated node pressure')
ylabel('Normalised SSE')
cb=colorbar
ylabel(cb,'Inside node fraction')

%%

figure(11)
clf

subplot(121)
hold on
box on
grid on

scatter(inv{1}.sse_fit,inv{1}.normsse_best3,10,inv{1}.inside_node_fraction,'filled')

xlabel('Entropy of estimated node pressure')
ylabel('Normalised SSE')
cb=colorbar
ylabel(cb,'Inside node fraction')


subplot(122)
hold on
box on
grid on

scatter(inv{2}.sse_fit,inv{2}.normsse_best3,10,inv{2}.inside_node_fraction,'filled')

xlabel('Entropy of estimated node pressure')
ylabel('Normalised SSE')
cb=colorbar
ylabel(cb,'Inside node fraction')

%%

figure(12)
clf

subplot(121)
hold on
box on
grid on


scatter(inv{1}.true_p_entropy_value,inv{1}.r_smooth,10,inv{1}.inside_node_fraction,'filled')
xlabel('Entropy of true node pressure')
ylabel('r')
cb=colorbar;
ylabel(cb,'Inside node fraction')

subplot(122)
hold on
box on
grid on


scatter(inv{1}.rl_entropy_value,inv{1}.r_smooth,10,inv{1}.inside_node_fraction,'filled')
xlabel('Entropy of RL')
ylabel('r')
cb=colorbar;
ylabel(cb,'Inside node fraction')

%%
figure(12)
clf
subplot(211)
hist(inv{1}.r_detailed,50)
subplot(212)
hist(inv{1}.r_smooth,50)
