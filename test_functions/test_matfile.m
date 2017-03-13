[f,d] = uigetfile();

mFile = matfile(fullfile(d,f),'writable',false);

tic
details = whos(mFile)
toc

