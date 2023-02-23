% Test using list of names for properties of a struct
s = struct();
fields = {'a', 'b', 'c'};
for k=1:3
   s.(fields{k}) = k*2;
end