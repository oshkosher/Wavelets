function AA=readDataFile(filename)
dataFl = fopen(filename,"r");
[sss,ccc,err] = fscanf(dataFl,"%s\n","C");
[dim,cnt] = fread(dataFl,2,"int32");
[ddd,cnt]=fread(dataFl,Inf,"float32");
data=reshape(ddd,dim);
fclose(dataFl);
figure;imagesc(data');colormap (gray(256));title(filename);axis equal;axis([1 dim(2) 1 dim(1)]);
AA = data';