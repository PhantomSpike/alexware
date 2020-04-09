function processed_data = cras(data)
%% Process the data to subtract the median
%Define the sourrounding channels which are used to compute and subtract
%the median
ix_mat{1} = [1:6];
ix_mat{2} = [1:7];
ix_mat{3} = [1:8];
ix_mat{4} = [1:9];
ix_mat{5} = [1:10];
ix_mat{380} = [375:384];
ix_mat{381} = [376:384];
ix_mat{382} = [377:384];
ix_mat{383} = [378:384];
ix_mat{384} = [379:384];

for ch_no = 6:379
    ix_mat{ch_no} = [ch_no - 5:ch_no + 5];
end

fprintf('== CRA processing ==\n');tic;
[rows,columns] = size(data);
median_data = data - median(data,2);
processed_data = zeros(rows,columns);
parfor ch_no = 1:384
    fprintf('== Channel no %.0f ==\n',ch_no);
    processed_data(ch_no,:) = median_data(ch_no,:) - median(median_data(ix_mat{ch_no},:));
end
processed_data(385,:) = data(385,:);
fprintf('== Done! Processing took %f sec ==\n',toc);
end