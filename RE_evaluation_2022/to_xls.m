ekf_data = load('grading/ekf_data.mat', 'data');
pf_data = load('grading/pf_data.mat', 'data');

writecell(ekf_data.data, 'grading/ekf_results.xls', 'UseExcel', false);
writecell(pf_data.data, 'grading/pf_results.xls', 'UseExcel', false);


%ekf_data = load('ekf_data.mat', 'data');
%writecell(ekf_data.data, 'ekf_results.xls', 'UseExcel', false);