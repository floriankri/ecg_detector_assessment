ver = 'v007';

root_data_folder = 'C:\Users\flori\OneDrive\Dokumente\TU\Bachelor Thesis\Code\Annotator\data\';
annotations_folder = [root_data_folder, '2022_annotations', filesep];
%csv_folder = [root_data_folder, 'csv', filesep];

database_path = 'C:\Users\flori\OneDrive\Dokumente\TU\Bachelor Thesis\Code\TestDatabases\annotator_test\';
% if ~exist(csv_folder, 'dir')
%         mkdir(csv_folder)
% end

files = dir(fullfile(annotations_folder,[ver,'*.mat']));

currentFolder = pwd;
cd(database_path)

for i = 1 : length(files)
    data = load([files(i).folder,filesep,files(i).name]);
    annotations = floor(data.pk_anns.t*data.up.fs);
    wrann([data.up.name],'atr',annotations);
%     if isempty(array)
%         array = -1;
%     end
%     dlmwrite([csv_folder,'test.csv'],transpose(array),'delimiter',',','-append');
end

cd(currentFolder);