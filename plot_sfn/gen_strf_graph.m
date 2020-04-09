function gen_strf_graph(animal,pen,neuron)
% animal = 'Ronnie';
% pen = 'P08';
% neuron = 11;
dir_name = ['/home/alex/Desktop/Kernels/',animal,'/',pen,'/data/sep_kernel_adjbeta.mat'];
save_dir = '/home/alex/Dropbox/SfN/';
coch_dir = '/mnt/40086D4C086D41D0/Reverb_data/coch.mat';
load(dir_name);
load(coch_dir);
%Generate STRF
anech_strf = sep_kernel(neuron).anech_k_h.k_f*sep_kernel(neuron).anech_k_h.k_h';
r1_strf = sep_kernel(neuron).reverb1_k_h.k_f*sep_kernel(neuron).reverb1_k_h.k_h';
r2_strf = sep_kernel(neuron).reverb2_k_h.k_f*sep_kernel(neuron).reverb2_k_h.k_h';
all_data = [anech_strf(:);r1_strf(:);r2_strf(:)];
mao = max(abs(all_data(:)));
anech_strf_n = anech_strf./mao;
r1_strf_n = r1_strf./mao;
r2_strf_n = r2_strf./mao;

%Plotting params
spacing_freq = 2;
spacing_time = 5;
freqs = coch(1).params.f;
freqs = fliplr(freqs);
num_freqs = numel(freqs);
freqs = ceil(freqs)/1000; %Convert the frequencies into kHz
freqs = freqs([1:spacing_freq:num_freqs]);
n_h = coch(1).params.n_h;
time_bin_ms = coch(1).params.time_bin_ms;
num_columns = 3;
num_rows = 1;
per = 0.02;
edgel = 0.06; edger = per; edgeh = per; edgeb = 0.07; space_h = 0.01; space_v =0.01;
[pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);
for ii = 1:numel(freqs)
    y_labels{ii} = num2str(freqs(ii),'%.1f');
end
time_ms = [-(n_h*time_bin_ms - time_bin_ms):spacing_time*time_bin_ms:0];

for jj = 1:numel(time_ms)
    x_labels{jj} = num2str(time_ms(jj),'%.0f');
end

%Actual plotting
%Anech
font1 = 20;
font2 = 20;
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('position',pos{1});
imagesc(anech_strf_n,[-1 1]);
colormap('redblue');
yticks([1:spacing_freq:num_freqs]);
yticklabels(y_labels);
set(gca,'FontName','Arial','FontSize',font1,'FontWeight','Bold');
ylabel('Frequency [kHz]','FontSize',font2,'FontWeight','bold');
xticks([1:spacing_time:n_h]);
xticklabels(x_labels);
xlabel('History [ms]','FontSize',font1,'FontWeight','bold');
% colorbar;
%Reverb 1
subplot('position',pos{2});
imagesc(r1_strf_n,[-1 1]);
colormap('redblue');
set(gca,'ytick',[]);
set(gca,'FontName','Arial','FontSize',font1,'FontWeight','Bold');
xticks([1:spacing_time:n_h]);
xticklabels(x_labels);
xlabel('History [ms]','FontSize',font1,'FontWeight','bold');
% colorbar;
%Reverb 2
subplot('position',pos{3});
imagesc(r2_strf_n,[-1 1]);
colormap('redblue');
set(gca,'ytick',[]);
set(gca,'FontName','Arial','FontSize',font1,'FontWeight','Bold');
xticks([1:spacing_time:n_h]);
xticklabels(x_labels);
xlabel('History [ms]','FontSize',font1,'FontWeight','bold');
colorbar;
save_name = [save_dir,animal,'_',pen,'_neuron',num2str(neuron),'.jpg'];
set(gcf,'color','w');
export_fig(save_name);
close;