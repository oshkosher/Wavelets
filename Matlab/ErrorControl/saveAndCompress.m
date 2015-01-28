function [ratio, original_bytes, compressed_bytes] = saveAndCompress(A, fileName)

% Setup global variables.
working_path = pwd;
COMPRESS_DIR = 'compress_dir';

% Generate the matrix file.
FILE=fileName;
%FMT='int';
[fid,mssg]=fopen(FILE,'w','l');
fwrite(fid, A);
status=fclose(fid);

% Get the matrix file size.
FILE_BYTES=strcat(working_path,'\',FILE);
s=dir(FILE_BYTES);
original_bytes = s.bytes;

% Compress the matrix file.
gzip(FILE,COMPRESS_DIR);

% Get the compressed original matrix file size.
FILE_BYTES=strcat(working_path,'\',COMPRESS_DIR,'\',FILE,'.gz');
s=dir(FILE_BYTES);
compressed_bytes = s.bytes;

ratio = original_bytes/compressed_bytes;