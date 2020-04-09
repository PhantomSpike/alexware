%% Plot the STRF
lag_max = 10;
resultsplot = results1x;
for n = 1:30
    results{n} = fliplr(resultsplot(n).kernel.k_fh);
    results{n} = flipud(results{n});
    results{n} = results{n}(:,1:lag_max);
end
freq = logspace(log10(170),log10(22000),35);
freq = fliplr(freq(1:30));
L = 50*(0:1:lag_max);

figure;
for n = 1:30
    subplot(6,5,n);
    maxabsk = max(max(abs(results{n}(:))));
    imagesc(L,freq,results{n},[-maxabsk maxabsk]);
    colormap(jet(100000));
    axis xy;
    caxis([-0.035 0.035]);
end
title=text(0.5,0.9,'Transformation kernels for room size 1x','FontSize',17,'FontWeight','bold');
xlab=text(0.5,0.1,'Time into the past [ms]','FontSize',17);
ylab=text(0.05,0.5,'Reverberant cochleagram frequency (low -> high) [Hz]','Rotation',90,'FontSize',17);
colorbar;
%% Ploting the model
model=model2x;
vec_z = zeros(1,25);
figure;
tbin = 50;
for ker_no = 1:30
    subplot(6,5,ker_no);
    plot(model{ker_no}(1,:)*tbin,model{ker_no}(2,:),'k.','markersize', 14);
    axis tight;
    hold on
    plot(model{ker_no}(1,:)*tbin,vec_z,'b--');
    plot(model{ker_no}(1,:)*tbin,model{ker_no}(3,:),'r','LineWidth',2);
    axis tight;
end
hold off;
title=text(0.1,0.1,'Time constants for room size 2x','FontSize',17,'FontWeight','bold');
xlab=text(0.5,0.1,'Time into the past [ms]','FontSize',17);
ylab=text(0.1,0.1,'Regression coefficient','Rotation',90,'FontSize',17);
%% Comparing the different rooms
freq = logspace(log10(170),log10(22000),35);
freq = freq(1:30);
figure;
semilogx(freq, betas1x,'LineWidth',4);
xlim([100 11000]);
hold on;
semilogx(freq, betas2x,'LineWidth',4);
xlim([100 11000]);
semilogx(freq, betas4x,'LineWidth',4);
xlim([100 11000]);
title('Frequency dependence of the time constants for rooms of different size');
xlabel('Frequency [Hz]');
ylabel('Time constant [ms]');
legend('Room size 1x','Room size 2x','Room size 4x');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',17);
hold off;

%% Plot the adaptation
freq = logspace(log10(170),log10(22000),35);
freq = freq(1:30);
figure;
semilogx(freq,1 - adapt1x,'LineWidth',4);
xlim([100 11000]);
hold on;
semilogx(freq,1 - adapt2x,'LineWidth',4);
xlim([100 11000]);
semilogx(freq, 1 - adapt4x,'LineWidth',4);
xlim([100 11000]);
title('Frequency dependence of the inhibition for rooms of different size');
xlabel('Frequency [Hz]');
ylabel('Amount of inhibition [% reduction of total]');
legend('Room size 1x','Room size 2x','Room size 4x');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',17);
hold off;