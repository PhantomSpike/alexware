fname = '/media/alex/5FC39EAD5A6AA312/GAD_labelling/analysis/';
save_dir = '/media/alex/5FC39EAD5A6AA312/GAD_labelling/Plots/';
%Parameters
tmt_mrk = 1; gad_mrk = 2; gcamp_mrk = 3; %These are the colour channel numbers that correspond to the different labels
epsilon = 7; %The minimum Euclidean distance for two labels to be considered the same
num_zplanes = 4;

dir_info = dir(fname);
dir_info = dir_info(~ismember({dir_info.name},{'.','..'}));


for z = 1:num_zplanes
    anna.data{z} = load(fullfile(dir_info(z).folder,dir_info(z).name)); %Load the data
    nala.data{z} = load(fullfile(dir_info(z+num_zplanes).folder,dir_info(z+num_zplanes).name));
    
    anna.gcamp_ix{z} = find(anna.data{z}.rez(:,7)==gcamp_mrk); %Find the ix of the coordinates that correspond to GCAMP labelled cells
    anna.gad_ix{z} = find(anna.data{z}.rez(:,7)==gad_mrk); %Find the ix of the coordinates that correspond to GAD labelled cells
    anna.tmt_ix{z} = find(anna.data{z}.rez(:,7)==tmt_mrk); %Find the ix of the coordinates that correspond to dTomato labelled cells
    
    anna.gcamp_xy{z} = anna.data{z}.rez(anna.gcamp_ix{z},5:6); %Extract the vectors corresponding to the label locations
    anna.gad_xy{z} = anna.data{z}.rez(anna.gad_ix{z},5:6);
    anna.tmt_xy{z} = anna.data{z}.rez(anna.tmt_ix{z},5:6);
    
    anna.gcamp_total(z) = length(anna.gcamp_xy{z});
    anna.gad_total(z) = length(anna.gad_xy{z});
    anna.tmt_total(z) = length(anna.tmt_xy{z});
    
    anna.label_mat{z} = zeros(anna.gcamp_total(z),3); %Initialize the label matrix where to keep track of the different labels
    anna.label_mat{z}(:,1) = 1; %Fill the GCAMP column with ones since we only care about cells that have GCAMP
    
    labels{1} = 'GCaMP';
    labels{2} = 'GAD';
    labels{3} = 'tdTomato';
 
    for n = 1:anna.gcamp_total(z)
        d_gad = anna.gad_xy{z} - anna.gcamp_xy{z}(n,:); %Find the difference vectors
        d_tmt = anna.tmt_xy{z} - anna.gcamp_xy{z}(n,:);
        
        dist_gad = sqrt(sum(d_gad.^2,2)); %Find the absolute difference between the two labels
        min_gad = min(dist_gad); %Find the smallest difference
        
        dist_tmt = sqrt(sum(d_tmt.^2,2));
        min_tmt = min(dist_tmt);
        
        anna.label_mat{z}(n,2) = min_gad<epsilon; %Fill in a 0 or 1 depending on whether the smallest distance is sufficinetly close
        anna.label_mat{z}(n,3) = min_tmt<epsilon;
    end
    
    categories{1} = 'tdTomato only';
    categories{2} = 'GAD only';
    categories{3} = 'tdTomato + GAD';
    categories{4} = 'GCaMP only';
    
    anna.labels_area(z,1) = sum(anna.label_mat{z}(:,2)==0 & anna.label_mat{z}(:,3)==1);
    anna.labels_area(z,2) = sum(anna.label_mat{z}(:,2)==1 & anna.label_mat{z}(:,3)==0);
    anna.labels_area(z,3) = sum(anna.label_mat{z}(:,2)==1 & anna.label_mat{z}(:,3)==1);
    anna.labels_area(z,4) = sum(anna.label_mat{z}(:,2)==0 & anna.label_mat{z}(:,3)==0);
    
    anna.labels_area_per(z,:) = 100*anna.labels_area(z,:)/sum(anna.labels_area(z,:)); %Convert to percentage 
    
    
    nala.gcamp_ix{z} = find(nala.data{z}.rez(:,7)==gcamp_mrk); %Find the ix of the coordinates that correspond to GCAMP labelled cells
    nala.gad_ix{z} = find(nala.data{z}.rez(:,7)==gad_mrk); %Find the ix of the coordinates that correspond to GAD labelled cells
    nala.tmt_ix{z} = find(nala.data{z}.rez(:,7)==tmt_mrk); %Find the ix of the coordinates that correspond to dTomato labelled cells
    
    nala.gcamp_xy{z} = nala.data{z}.rez(nala.gcamp_ix{z},5:6); %Extract the vectors corresponding to the label locations
    nala.gad_xy{z} = nala.data{z}.rez(nala.gad_ix{z},5:6);
    nala.tmt_xy{z} = nala.data{z}.rez(nala.tmt_ix{z},5:6);
    
    nala.gcamp_total(z) = length(nala.gcamp_xy{z});
    nala.gad_total(z) = length(nala.gad_xy{z});
    nala.tmt_total(z) = length(nala.tmt_xy{z});
    
    nala.label_mat{z} = zeros(nala.gcamp_total(z),3); %Initialize the label matrix where to keep track of the different labels
    nala.label_mat{z}(:,1) = 1; %Fill the GCAMP column with ones since we only care about cells that have GCAMP
    
    for n = 1:nala.gcamp_total(z)
        d_gad = nala.gad_xy{z} - nala.gcamp_xy{z}(n,:); %Find the difference vectors
        d_tmt = nala.tmt_xy{z} - nala.gcamp_xy{z}(n,:);
        
        dist_gad = sqrt(sum(d_gad.^2,2)); %Find the absolute difference between the two labels
        min_gad = min(dist_gad); %Find the smallest difference
        
        dist_tmt = sqrt(sum(d_tmt.^2,2));
        min_tmt = min(dist_tmt);
        
        nala.label_mat{z}(n,2) = min_gad<epsilon; %Fill in a 0 or 1 depending on whether the smallest distance is sufficinetly close
        nala.label_mat{z}(n,3) = min_tmt<epsilon;
    end
    
    nala.labels_area(z,1) = sum(nala.label_mat{z}(:,2)==0 & nala.label_mat{z}(:,3)==1);
    nala.labels_area(z,2) = sum(nala.label_mat{z}(:,2)==1 & nala.label_mat{z}(:,3)==0);
    nala.labels_area(z,3) = sum(nala.label_mat{z}(:,2)==1 & nala.label_mat{z}(:,3)==1);
    nala.labels_area(z,4) = sum(nala.label_mat{z}(:,2)==0 & nala.label_mat{z}(:,3)==0);
    
    nala.labels_area_per(z,:) = 100*nala.labels_area(z,:)/sum(nala.labels_area(z,:)); %Convert to percentage 
end

%% Plotting
%First plot the individual areas per animal
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w');
precision = 3;
row = 2;
col = 4;
per = 0.04;
edgel = 0.1; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = 0.01;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
for jj = 1:num_zplanes
    subplot('position',pos{jj});
    labels = {['Tomato',newline,'+GFP ',newline,num2str(anna.labels_area_per(jj,1),precision),'%'],...
        ['GAD',newline,'+GFP ',num2str(anna.labels_area_per(jj,2),precision),'%'],...
        ['Tomato',newline,'+GAD',newline,'+GFP ',newline,num2str(anna.labels_area_per(jj,3),precision),'%'],...
        ['GFP ',newline,num2str(anna.labels_area_per(jj,4),precision),'%']};
    p = pie(anna.labels_area(jj,:),labels);
    title(['Anna Area ',num2str(jj),' n=',num2str(anna.gcamp_total(jj))]);
    set(gca,'FontWeight','bold')
    colormap([1 0 0;      %// red
        0 0 1;      %// green
        1 0 1;      %// blue
        0 1 0])  %// grey
    
    subplot('position',pos{num_zplanes+jj});
    labels = {['Tomato',newline,'+GFP',newline,num2str(nala.labels_area_per(jj,1),precision),'%'],...
        ['GAD',newline,'+GFP',newline,num2str(nala.labels_area_per(jj,2),precision),'%'],...
        ['Tomato',newline,'+GAD',newline,'+GFP',newline,num2str(nala.labels_area_per(jj,3),precision),'%'],...
        ['GFP',newline,num2str(nala.labels_area_per(jj,4),precision),'%']};
    p = pie(nala.labels_area(jj,:),labels);
    title(['Nala Area ',num2str(jj),' n=',num2str(nala.gcamp_total(jj))]);
    colormap([1 0 0;      %// red
        0 0 1;      %// green
        1 0 1;      %// blue
        0 1 0])  %// grey
end
save_name = [save_dir,'Cell_quant_Anna_Nala_areas.png'];
export_fig(save_name);
close all;

%Plot the average from all areas together
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w');
precision = 3;
row = 1;
col = 2;
per = 0.04;
edgel = 0.1; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = 0.01;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
subplot('position',pos{1});
anna_mean = mean(anna.labels_area_per);
labels = {['Tomato',newline,'+GFP ',newline,num2str(anna_mean(1),precision),'%'],...
    ['GAD',newline,'+GFP ',num2str(anna_mean(2),precision),'%'],...
    ['Tomato',newline,'+GAD',newline,'+GFP ',newline,num2str(anna_mean(3),precision),'%'],...
    ['GFP ',newline,num2str(anna_mean(4),precision),'%']};
p = pie(anna_mean,labels);
title(['Anna Total ',' n=',num2str(sum(anna.gcamp_total))]);
set(gca,'FontWeight','bold')
colormap([1 0 0;      %// red
    0 0 1;      %// green
    1 0 1;      %// blue
    0 1 0])  %// grey

subplot('position',pos{2});
nala_mean = mean(nala.labels_area_per);
labels = {['Tomato',newline,'+GFP ',newline,num2str(nala_mean(1),precision),'%'],...
    ['GAD',newline,'+GFP ',num2str(nala_mean(2),precision),'%'],...
    ['Tomato',newline,'+GAD',newline,'+GFP ',newline,num2str(nala_mean(3),precision),'%'],...
    ['GFP ',newline,num2str(nala_mean(4),precision),'%']};
p = pie(nala_mean,labels);
title(['Nala Total ',' n=',num2str(sum(nala.gcamp_total))]);
set(gca,'FontWeight','bold')
colormap([1 0 0;      %// red
    0 0 1;      %// green
    1 0 1;      %// blue
    0 1 0])  %// grey
save_name2 = [save_dir,'Cell_quant_Anna_Nala_Total.png'];
export_fig(save_name2);
close all;