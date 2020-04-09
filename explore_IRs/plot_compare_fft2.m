function plot_compare_fft2(data1,data2,meta)

font_axis = 12;
line_width = 2;
fs = meta.fs;

[P1,f] = fft_run(data1,fs,false);
[P2,f] = fft_run(data2,fs,false);
full = 'full';
if meta.short_fft ~= 0
    full = 'short';
    [~, ix_match] = min(abs(f-meta.short_fft/1000));
    window_ix = meta.f_window./(mean(diff(f))*1000); %Find the ix for the frqeuncy plotting window
    f = f(ix_match-window_ix:ix_match+window_ix);
    P1 = P1(ix_match-window_ix:ix_match+window_ix);
    P2 = P2(ix_match-window_ix:ix_match+window_ix);
end

max_P = max([P1(:);P2(:)]);
figure('units','normalized','outerposition',[0 0 1 1]);
row = 1;
col = 2;
per = 0.005;
edgel = 0.03; edger = per; edgeh = 0.03; edgeb = 0.05; space_h = per; space_v = per;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
subplot('position',pos{1});
plot(f,P1,'LineWidth',line_width);
ylim([0 max_P]);
legend([meta.var_name,' normal ',meta.ear,' ear']);
xlabel('Freqeuncy [kHz]');
ylabel('Fourier magnitude');
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
subplot('position',pos{2});
plot(f,P2,'Color',[0.8500 0.3250 0.0980],'LineWidth',line_width);
ylim([0 max_P]);
set(gca, 'YTick', []);
legend([meta.var_name,' bp between ',meta.lowf,'-',meta.highf,'Hz ',meta.ear,' ear']);
xlabel('Freqeuncy [kHz]');
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
set(gcf,'color','w');
save_name = [meta.save_dir,meta.var_name,'_',full,'_normal_vs_bp_',meta.lowf,'-',meta.highf,'Hz_','_FFT_spectrum_',meta.ear,'_ear.jpg'];
export_fig(save_name);
close all;