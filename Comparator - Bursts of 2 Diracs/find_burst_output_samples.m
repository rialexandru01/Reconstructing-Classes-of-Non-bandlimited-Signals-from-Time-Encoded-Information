function [ t_stop_burst_estimation, y_n_burst, n_vec_burst, y_n_temp ] = find_burst_output_samples( y_n_temp, t_stop_burst_estimation, y_n1, n_vec1 )
%FIND_BURST_OUTPUT_SAMPLES Find the output samples corresponding to the
%current burst

t_start_burst_estimation = find(y_n_temp ~= 0, 1, 'first') + t_stop_burst_estimation;
y_n_temp = y_n_temp(find(y_n_temp ~= 0, 1, 'first'):end);
t_stop_burst_estimation = t_start_burst_estimation + find(y_n_temp == 0, 1, 'first')-2;
y_n_temp = y_n_temp(find(y_n_temp == 0, 1, 'first'):end);
y_n_burst = y_n1(t_start_burst_estimation:t_stop_burst_estimation);
n_vec_burst = n_vec1(t_start_burst_estimation:t_stop_burst_estimation);


end

