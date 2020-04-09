animal = 'Ronnie';
pen = 'P13';
neuron = 7;
dir_name = ['/home/alex/Desktop/Kernels/',animal,'/',pen,'/data/sep_kernel_adjbeta.mat'];
save_dir = '/home/alex/Dropbox/SfN/';
load(dir_name);
spacing_time = 5;
n_h = 30;
time_bin_ms = 10;
time_ms = [-(n_h*time_bin_ms - time_bin_ms):time_bin_ms:0];
% time_ms = [-(n_h*time_bin_ms - time_bin_ms):spacing_time*time_bin_ms:0];
% for jj = 1:numel(time_ms)
%     x_labels{jj} = num2str(time_ms(jj),'%.0f');
% end

k_h.anech = sep_kernel(neuron).anech_k_h.k_h; 
k_h.r1 = sep_kernel(neuron).reverb1_k_h.k_h; 
k_h.r2 = sep_kernel(neuron).reverb2_k_h.k_h; 
mao = max(abs([k_h.anech;k_h.r1;k_h.r2]));
k_h.anech = k_h.anech./mao; 
k_h.r1 = k_h.r1./mao;
k_h.r2 = k_h.r2./mao;
F.anech = sep_kernel(neuron).dexp_anech.F;
F.r1 = sep_kernel(neuron).dexp_reverb1.F;
F.r2 = sep_kernel(neuron).dexp_reverb2.F;
F.anech = F.anech./mao;
F.r1 = F.r1./mao;
F.r2 = F.r2./mao;


%Plotting
sz = 30;
line_sz = 5;
font_y = 20;
font_x = 20;
font_axis = 16;
num_columns = 3;
num_rows = 1;
per = 0.02;
edgel = 0.07; edger = per; edgeh = 0.05; edgeb = 0.07; space_h = 0.03; space_v =0.01;
[pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('position',pos{1});
hold on;
plot(time_ms,k_h.anech,'k.','MarkerSize',sz);
time_ms_plot = time_ms(1+numel(k_h.anech) - numel(F.anech):end);
plot(time_ms_plot,fliplr(F.anech));
ylim([-1 1]);
hold off;
% axis tight;
xlabel('History [ms]','FontSize',font_x,'FontWeight','bold');
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
ylabel('Weight value [AU]','FontSize',font_y,'FontWeight','bold');
set(findall(gca, 'Type', 'Line'),'LineWidth',line_sz);

subplot('position',pos{2});
hold on;
plot(time_ms,k_h.r1,'k.','MarkerSize',sz);
time_ms_plot = time_ms(1+numel(k_h.r1) - numel(F.r1):end);
plot(time_ms_plot,fliplr(F.r1),'k');
ylim([-1 1]);
set(gca,'ytick',[]);
hold off;
% axis tight;
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
xlabel('History [ms]','FontSize',font_x,'FontWeight','bold');
set(findall(gca, 'Type', 'Line'),'LineWidth',line_sz);

subplot('position',pos{3});
hold on;
plot(time_ms,k_h.r2,'k.','MarkerSize',sz);
time_ms_plot = time_ms(1+numel(k_h.r2) - numel(F.r2):end);
plot(time_ms_plot,fliplr(F.r2),'r');
ylim([-1 1]);
set(gca,'ytick',[]);
hold off;
% axis tight;
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
xlabel('History [ms]','FontSize',font_x,'FontWeight','bold');
set(findall(gca, 'Type', 'Line'),'LineWidth',line_sz);
set(gcf,'color','w');
save_name = [save_dir,animal,'_',pen,'_neuron',num2str(neuron),'_dexp_.jpg'];
export_fig(save_name);
close;