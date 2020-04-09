%Script to get the mean lfps when noise bursts are presented
close all
clear all

rootdir='F:\Ioanna\22_03_2019\P1\Gaps';
penetration = rootdir(23);

lfpD = dir(fullfile(rootdir,  '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(rootdir , lfpD(1).name);
load([rootdir '\meta_info\synch_ch.mat'])
fs = 30000; %sampling rate of the probe
lfpFs = 2500;  % neuropixels phase3a
min_trig_length_s = 0.5; %sec
min_inter_trig_length_s= 0.5; %sec
[start_time_ms] = get_triggers_new(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);
start_time_lfpFs=round(start_time_ms./1000.*lfpFs); %convert stim start times to lfp samples

nChansInFile = 385;  % neuropixels phase3a, from spikeGLX. 385 or 277

%load lfp data
lfpdata=getLFPdata(lfpFilename, nChansInFile);

resp_dur=150; %dur of clip to take, in ms
resp_dur_lfpFs=resp_dur./1000.*lfpFs; %convert times to lfp samples
offset_dur=30; %start lfp clip by this amount before stimulus, in ms
offset_dur_lfpFs=offset_dur./1000.*lfpFs; %convert times to lfp samples
lfp_responses=[];
mean_lfp_responses=zeros(nChansInFile,resp_dur_lfpFs);

for ii=1:length(start_time_lfpFs),
    lfp_responses(ii).traces=lfpdata(:,(start_time_lfpFs(ii)-offset_dur_lfpFs):(start_time_lfpFs(ii)-offset_dur_lfpFs+resp_dur_lfpFs-1));
    mean_lfp_responses=mean_lfp_responses+lfp_responses(ii).traces;
end
mean_lfp_responses=mean_lfp_responses./length(start_time_lfpFs);
mean_lfp_responses=flipud(mean_lfp_responses);

figure(1);clf
imagesc(mean_lfp_responses)
colormap(flipud(jet))
caxis([-40 15])


xt=get(gca,'xtick');
set(gca,'xticklabel',xt./lfpFs.*1000)
hx = xlabel('Time (ms)');
yt=get(gca,'ytick');
set(gca,'yticklabel',yt*10)
hy = ylabel('Depth on the probe (µm)');
ht = title(strcat('Penetration',penetration,' mean LFP'));
set(gca,'FontName','Georgia','FontSize',8)

set([ht,hy,hx],'FontName','Georgia')
set([hy,hx],'FontSize',14)
set(ht,'FontSize',14)

colorbar('Ticks' ,[-40 0 15],'TickLabels' ,{ '-', '0', '+'},'FontName','Georgia','FontSize',12)

p=[1295         229         634         767];
set(gcf,'position',p);


%% script to add the lines and calculate the channel and the location limits
%of our layers


%firtst set the parameters:
starting_point = 0; %location of brain surface?123
reversal_point=100 ; %?m of location of lfp reversal point
l23_l4_limit=300; %?m of location of limit between L2/3 & L4
l4_l5_limit=500; %?m of location of limit between L4 & L5
cortical_limit=1000; %?m of location of lower cortical limit
chan_dist=10; %assumed distance between channels to calculate
%the channel to distance transformation

%then add the y location of the image where you set your reversal point manually

upper_point = 213;


hold on

%then calculate the limits coordinates
upper_cortical_limit_y_coord = upper_point;
%coordinate equivalent of upper cortical limit

reversal_point_coord = upper_point + reversal_point/chan_dist;

limit_l23_l4_y_coord = upper_point + l23_l4_limit/chan_dist;
%coordinate equivalent of limit between superficial and middle layers


limit_l4_l5_y_coord = upper_point + l4_l5_limit/chan_dist;
%coordinate equivalent of limit between middle and deep layers

lower_cortical_limit_y_coord = upper_point + cortical_limit/chan_dist;
%coordinate equivalent of lower cortical limit

hold on

line([0 size(mean_lfp_responses,2)],[reversal_point_coord reversal_point_coord],'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')

line([0 size(mean_lfp_responses,2)],[upper_cortical_limit_y_coord upper_cortical_limit_y_coord ],'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
line([0 size(mean_lfp_responses,2)],[limit_l23_l4_y_coord  limit_l23_l4_y_coord ],'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
line([0 size(mean_lfp_responses,2)],[limit_l4_l5_y_coord limit_l4_l5_y_coord],'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
line([0 size(mean_lfp_responses,2)],[lower_cortical_limit_y_coord lower_cortical_limit_y_coord],'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')

text(size(mean_lfp_responses,2)/25 ,limit_l4_l5_y_coord+(lower_cortical_limit_y_coord-limit_l4_l5_y_coord)/2, 'Layer 5-6','Color','k','FontSize',12,'FontName','Georgia')
text(size(mean_lfp_responses,2)/25 ,limit_l23_l4_y_coord+(limit_l4_l5_y_coord-limit_l23_l4_y_coord)/2, 'Layer 4','Color','k','FontSize',12,'FontName','Georgia')
text(size(mean_lfp_responses,2)/25 ,reversal_point_coord+(limit_l23_l4_y_coord-reversal_point_coord)/2, 'Layer 2-3','Color','k','FontSize',12,'FontName','Georgia')
text(size(mean_lfp_responses,2)/25 ,upper_cortical_limit_y_coord+(reversal_point_coord-upper_cortical_limit_y_coord)/2, 'Layer 1','Color','k','FontSize',10,'FontName','Georgia')

limits=[upper_cortical_limit_y_coord,limit_l23_l4_y_coord,limit_l4_l5_y_coord,lower_cortical_limit_y_coord];

channel_limits = zeros(size(limits));

for ii = 1:length(limits)
    channel_limits(ii) = nChansInFile - limits(ii)+1;
end

limits_depths = channel_limits.*chan_dist;

save(['C:\Users\Ioanna\Desktop\10_04_2019\Results\',num2str(penetration),'depth_limits.mat'],'limits_depths')
save(['C:\Users\Ioanna\Desktop\10_04_2019\Results\',num2str(penetration),'upper_channel.mat'],'upper_point')

disp(['Layers 1/2/3 range from  ' num2str(limits_depths(1)) ' till ' num2str(limits_depths(2)) ' um on the probe'])

disp(['Layer 4 ranges from  ' num2str(limits_depths(2)) ' till ' num2str(limits_depths(3)) ' um on the probe' ])

disp(['Layers 5/6 range from  ' num2str(limits_depths(3)) ' till ' num2str(limits_depths(4)) ' um on the probe' ])

hold on
load F:\Ioanna\22_03_2019\P1\Clicks\clust_info.mat %%load cluster info
load('C:\Users\Ioanna\Desktop\10_04_2019\Results\significant_all_rates_P1.mat')
plot_position = (384- round(clust_info.clust_depth./10)+1);
hp = plot(zeros(size(plot_position))+length(mean_lfp_responses)-60,plot_position,'o');
set(hp,'Marker', 'o', 'MarkerSize', 4,'MarkerEdgeColor' , 'k','MarkerFaceColor' , 'k');

plot_position = (384- round(clust_info.clust_depth(significant)./10)+1);
hp = plot(zeros(size(plot_position))+length(mean_lfp_responses)-50,plot_position,'o');
set(hp,'Marker', 'o', 'MarkerSize', 4,'MarkerEdgeColor' , 'r','MarkerFaceColor' , 'r');

save_dir = 'C:\Users\Ioanna\Desktop\10_04_2019\Results\';
print(strcat(save_dir, num2str(penetration),'lfps', '.png'),'-dpng')
savefig([save_dir, num2str(penetration),'lfps','.fig'])



%%
close all
%plot traces around region of proposed reversal point

figure(1);clf
plot_ch=[upper_point-5:4:upper_point+120];
for ii=1:length(plot_ch),
    plot(mean_lfp_responses(plot_ch(ii),:)-ones(1,size(mean_lfp_responses(plot_ch(ii),:),2)).*ii*10);
    hold on
end

hp = plot(mean_lfp_responses(upper_point,:)-ones(1,size(mean_lfp_responses(upper_point,:),2)).*2.5*10);
set(hp, 'LineWidth',2.5)
axis([0 350 -300 0])
xt=get(gca,'xtick');
set(gca,'xticklabel',xt./lfpFs.*1000)
hx = xlabel('Time (ms)');
yt=get(gca,'ytick');
set(gca,'yticklabel',[ 3240        3040        2840        2640        2440        2240        2040])
hy = ylabel('Depth on the probe');
ht = title(strcat('Penetration',penetration));
set(gca,'FontName','Georgia','FontSize',8)
p=[1295         229         634         767];
set(gcf,'position',p);


set([ht,hy,hx],'FontName','Georgia')
set([hy,hx],'FontSize',12)
set(ht,'FontSize',14)

save_dir = 'C:\Users\Ioanna\Desktop\10_04_2019\Results\';
print(strcat(save_dir, num2str(penetration),'channels', '.png'),'-dpng')
savefig([save_dir, num2str(penetration),'channels','.fig'])

