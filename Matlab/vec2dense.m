function vec2dense( v, filename )

str=strcat(filename,'.txt');
fid=fopen(str,'w');if ( fid==-1 ) error('could not open file'); end;
fprintf(fid,'%u\n',length(v)');
if(isinteger(v))
    fprintf(fid,'%d\n',v);
else
    fprintf(fid,'%.16f\n',v);
end
fclose(fid);