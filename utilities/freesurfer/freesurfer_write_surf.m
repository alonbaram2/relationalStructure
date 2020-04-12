function [vertex_coords, faces] = freesurfer_write_surf(fname,vertex,faces)
% function [vertex_coords, faces] = freesurfer_write_surf(fname,vertex,faces)
% Writes out a freesurfer surface file 
% vertex-indices for faces have to be 0-based! 
% Only accepts triagular meshes 
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
% 

vnum=size(vertex,1); 
fnum=size(faces,1); 

TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
QUAD_FILE_MAGIC_NUMBER =  16777215 ;

fid = fopen(fname, 'wb', 'b') ;
if (fid < 0)
  str = sprintf('Error opening file %s.', fname) ;
  error(str) ;
end
fwrite3(fid,TRIANGLE_FILE_MAGIC_NUMBER) ;
fprintf(fid,'New surface\n\n'); 
fwrite(fid,vnum,'int32') ;
fwrite(fid,fnum,'int32') ;
fwrite(fid,vertex','float32'); 
fwrite(fid, faces', 'int32') ;
fclose(fid) ;
