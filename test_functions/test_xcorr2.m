for n=1:1000
    
    c1 = conv2(rand(64,84),rand(64,128));
    
    c2 = xcorr2(rand(64,84),rand(64,128));
end