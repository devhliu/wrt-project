% Reading a function on a grid [ngrid x ngrid x ngrid] from file
% It is assumed that function is given in a CSV file with ordering 
% x -> y -> z (starting with a minimal index)

function grid_function = read_func_grid3d(ngrid, file)
  grid_function = zeros(ngrid, ngrid, ngrid);
  function_array = csvread(file);
  for i = 1 : ngrid
    for j = 1 : ngrid
      for k = 1 : ngrid
       grid_function(ngrid - j + 1, i, k) = init_point_data((i-1)*ngrid^2 + (j-1)*ngrid + k); 
      end
    end
  end
  % return grid_function
end 
