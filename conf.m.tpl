%% startup configuration
clear all
% exec_dir = '/usr/local/MATLAB/R2012a/mmder/plantjaja/';
exec_dir = '{exec_dir}';
% upload_dir = '/usr/local/MATLAB/R2012a/mmder/plantjaja/TestFolder';
upload_dir = '{upload_dir}';
% filename = 'IM001.jpg';
filename = '{filename}';
full_path = strcat(upload_dir, '/', filename);
% phyllotaxy = 1; % the definition/type of leaf look (single leaf or multiple leaves?)
phyllotaxy = {phyllotaxy}; % the definition/type of leaf look (single leaf or multiple leaves?)
