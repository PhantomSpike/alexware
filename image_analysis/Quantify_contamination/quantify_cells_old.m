filename = '/home/alex/Desktop/Code/alexware/image_analysis/Quantify_contamination/Nala1_4.mat';
load(filename);
tmt_mrk = 1; gad_mrk = 2; gcamp_mrk = 3; %These are the colour channel numbers that correspond to the different labels
epsilon = 3; %The minimum Euclidean distance for two labels to be considered the same

gcamp_ix = find(rez(:,7)==gcamp_mrk); %Find the ix of the coordinates that correspond to GCAMP labelled cells
gad_ix = find(rez(:,7)==gad_mrk); %Find the ix of the coordinates that correspond to GAD labelled cells
tmt_ix = find(rez(:,7)==tmt_mrk); %Find the ix of the coordinates that correspond to dTomato labelled cells

gcamp_xy = rez(gcamp_ix,5:6); %Extract the vectors corresponding to the label locations
gad_xy = rez(gad_ix,5:6);
tmt_xy = rez(tmt_ix,5:6);

num_neurons = length(gcamp_xy);
label_mat = zeros(num_neurons,3); %Initialize the label matrix where to kepe track fo the different labels
label_mat(:,1) = 1; %Fill the GCAMP column with ones since we only care about cells that have GCAMP

for n = 1:num_neurons
    d_gad = gad_xy - gcamp_xy(n,:); %Find the difference vectors
    d_tmt = tmt_xy - gcamp_xy(n,:);
    dist_gad = sqrt(sum(d_gad.*d_gad,2)); %Find the absolute difference between the two labels
    min_gad = min(dist_gad); %Find the smallest difference
    dist_tmt = sqrt(sum(d_tmt.*d_tmt,2));
    min_tmt = min(dist_tmt);
    
    label_mat(n,2) = min_gad<epsilon; %Fill in a 0 or 1 depending on whether the smallest distance is sufficinetly close
    label_mat(n,3) = min_tmt<epsilon;
end

marker_size = 30;
figure;
hold on;
title('Normal colors');
xlabel('X pos [um]');
ylabel('Y pos [um]');
plot(gcamp_xy(:,1),gcamp_xy(:,2),'g.','MarkerSize',marker_size);
plot(gad_xy(:,1),gad_xy(:,2),'b.','MarkerSize',marker_size);
plot(tmt_xy(:,1),tmt_xy(:,2),'r.','MarkerSize',marker_size);
hold off;

figure;
hold on;
title('Labels');
xlabel('X pos [um]');
ylabel('Y pos [um]');
both_ix = find(label_mat(:,2)==1 & label_mat(:,3)==1);
gad_only_ix = find(label_mat(:,2)==1 & label_mat(:,3)==0);
tmt_only_ix = find(label_mat(:,2)==0 & label_mat(:,3)==1);
no_lbl_ix = find(label_mat(:,2)==0 & label_mat(:,3)==0);

plot(gcamp_xy(both_ix,1),gcamp_xy(both_ix,2),'m.','MarkerSize',marker_size);
plot(gcamp_xy(gad_only_ix,1),gcamp_xy(gad_only_ix,2),'b.','MarkerSize',marker_size);
plot(gcamp_xy(tmt_only_ix,1),gcamp_xy(tmt_only_ix,2),'r.','MarkerSize',marker_size);
plot(gcamp_xy(no_lbl_ix,1),gcamp_xy(no_lbl_ix,2),'k.','MarkerSize',marker_size);
hold off;

    
    