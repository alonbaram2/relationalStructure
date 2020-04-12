function R=caret_getregions(varargin) 
% GUI interface for creating regions 
%           R.type='paint'
%           R.location=[filename];
%           --------------------
%           R.type='node'
%           R.location=nodenum
types={'done','paint','node'};
R={};p='x';r=1;
while (~strcmp(p,'done'))
    [p,YPos] = spm_input('Select region type',1,'m',types,types,1);
    type=p{1};
    switch (type)
        case 'paint'
            R{r}.type=type;
            R{r}.location=spm_get(inf,'*.paint','Select ROI images');
        case 'node'
            R{r}.type=type;
            R{r}.location=spm_input('Node Number','+1','e');
    end;
    r=r+1;
end;
    