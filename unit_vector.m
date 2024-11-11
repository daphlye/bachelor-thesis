function [output] = unit_vector(vector_array)
%UNIT_VECTOR calculates the direction of a vector
%   Input can be both an array with many vectors or just one vector
output = vector_array./vecnorm(vector_array, 2, 2);
end

