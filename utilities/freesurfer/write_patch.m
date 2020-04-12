%
% write_patch.m
%


function write_patch(fname,data,indx)


fid = fopen(fname,'w','b');
if (fid == -1)
   error('could not open file %s', fname) ;
end

fwrite(fid, -1,'int'); 

npts=size(indx,1); 
fwrite(fid, npts, 'int',0,'b') ;

for i=1:npts
    fwrite(fid, indx(i), 'int',0,'b') ;
    fwrite(fid, data(i,1), 'float',0,'b') ;
    fwrite(fid, data(i,2), 'float',0,'b') ;
    fwrite(fid, data(i,3), 'float',0,'b') ;
end

fclose(fid);
