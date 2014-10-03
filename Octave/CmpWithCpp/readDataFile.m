% XX = readDataFile('chicago_orig.data');
function AA=readDataFile(filename)
dataFl = fopen(filename,'r');
sss = fscanf(dataFl,'%s\n',1);
if ( strcmp(sss,'binary') == 0)
	display('Wrong type!  You need a binary .data file')
	AA = 0;
else
    [dim,cnt1] = fread(dataFl,2,'int32');
    dim = reshape(dim,[1 2]);
    [ddd,cnt2] = fread(dataFl,Inf,'float32');
    data=reshape(ddd,dim);
    fclose(dataFl);
    figure;imagesc(data');colormap (gray(256));title(filename);axis equal;axis([1 dim(2) 1 dim(1)]);
    AA = data';
    fdim = size(AA);
    if (fdim ~= dim)
        display('Dimentions of the array read don''t match the header!')
    end

end
