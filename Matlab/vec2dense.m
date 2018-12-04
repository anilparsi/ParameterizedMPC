function vec2dense( v, filename )
% Saves the contents of the vector v in the file filename.txt. The first 
% line is the number of elements in v. 

str=strcat(filename,'.txt');
fid=fopen(str,'w');

if ( fid==-1 ) 
    error('could not open file'); 
end

fprintf(fid,'%u\n',length(v)');

if(isinteger(v))
    fprintf(fid,'%d\n',v);
else
    fprintf(fid,'%.16f\n',v);
end

fclose(fid);
end