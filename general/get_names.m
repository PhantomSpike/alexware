function [file_name,folder_name] = get_names(file_path,rep)
%This function returns the filename of a given file from the full file_path
%, removing all the folder names before the actual file. It also replaces underscores '_' with empty spaces '
%'. Additionally, it also provides the folder name that the files is in.
%There is an extra option to remove underscores which is useful for names
%used in plotting
%>>INPUT>>
%file_path = The absolute path to a given file
%rep(optional) = A Boolean that determines whether to replaces underscores
%in the name with empty spaces. This is only useful for getting the names
%for plotting
%<<OUTPUT<<
%file_name = The name of the file without the extensions
%folder_name = The name of the folder containing the file

if nargin < 2
    rep = 0;
end

slash_ix = find(file_path == '/', 1, 'last');
dot_ix = find(file_path == '.', 1, 'last');
file_name = file_path(slash_ix+1:dot_ix-1);

if rep
    file_name = strrep(file_name,'_',' ');
end

folder_name = file_path(1:slash_ix);