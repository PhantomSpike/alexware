function plot_compare_coch(file_path1,file_path2,save_dir,ear)

time_step = 18;
freq_step = 3;

if strcmp(ear,'left')
    chan = 1;
elseif strcmp(ear,'right')
    chan = 2;
end

slash_ix1 = find(file_path1 == '/', 1, 'last');
slash_ix2 = find(file_path2 == '/', 1, 'last');
var_name1 = file_path1(slash_ix1+1:end);
var_name2 = file_path2(slash_ix2+1:end);
var_name1 = strrep(var_name1,'_',' ');
var_name2 = strrep(var_name2,'_',' ');


stim1 = load(file_path1);
stim2 = load(file_path2);
fs = stim1.Fs;

% [stim1.data,fs] = audioread(file_path1);
% [stim2.data,fs] = audioread(file_path2);


data1 = stim1.data(:,chan);
data2 = stim2.data(:,chan);
sz1 = length(data1);
sz2 = length(data2);

if sz1 > sz2
    data2(sz2:sz1) = 0;
else
    data1(sz1:sz2) = 0;
end

plot_on = false;
[X_ft1] = example_coch(data1,fs,plot_on);
[X_ft2,t,params] = example_coch(data2,fs,plot_on);

min_val = min([X_ft1(:);X_ft2(:)]);
max_val = max([X_ft1(:);X_ft2(:)]);

row = 1;
col = 2;
per = 0.005;
edgel = 0.05; edger = per; edgeh = 0.03; edgeb = 0.05; space_h = per; space_v = per;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('position',pos{1});
imagesc(X_ft1);
time_s = t(1:time_step:numel(t));
xticks(1:time_step:numel(t));
for jj = 1:numel(time_s)
    x_labels{jj} = num2str(time_s(jj),'%.2f');
end
xticklabels(x_labels);

yticks(1:freq_step:numel(params.f));
freqs = params.f(1:freq_step:numel(params.f))/1000; %Convert the frequencies into kHz
for jj = 1:numel(freqs)
    y_labels{jj} = num2str(freqs(jj),'%.2f');
end
yticklabels(fliplr(y_labels));
set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold');
xlabel('Time [s]','FontSize',12,'FontWeight','bold');
ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
colormap('inferno');
caxis([min_val max_val]);
subplot('position',pos{2});
imagesc(X_ft2);
xticks(1:time_step:numel(t));
xticklabels(x_labels);
yticklabels(fliplr(y_labels));
set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold');
xlabel('Time [s]','FontSize',12,'FontWeight','bold');
set(gca, 'YTick', []);
colormap('inferno');
% colorbar;
caxis([min_val max_val]);
set(gcf,'color','w');
save_name = [save_dir,var_name2,'_vs_',var_name1,'_cochleagram_',ear,'_ear.jpg'];
export_fig(save_name);
close all;
