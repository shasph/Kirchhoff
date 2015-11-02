function [V,F] = readobj()
  % Read the obj mesh file to save it to mat file
  
  V = zeros(0,3);
  F = zeros(0,3);
  vertex_index = 1;
  face_index = 1;
  fid = fopen('plate2.obj','rt');
  line = fgets(fid);
  while ischar(line)
    vertex = sscanf(line,'v %f %f %f');
    face = sscanf(line,'f %d %d %d');

    if(size(vertex)>0)
      V(vertex_index,:) = vertex;
      vertex_index = vertex_index+1;
    elseif(size(face,1)==3)
      F(face_index,:) = face;
      face_index = face_index+1;
     end
    save('plate2.mat','V','F');
    line = fgets(fid);   
  end
  fclose(fid);
end