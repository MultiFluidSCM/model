function field_out = initial_field( z_in, field_z, field )

n = length(field_z);
if n == 1
    field_out = field;
else
    i = 2;
    while z_in > field_z(i) && i < n
        i = i + 1;
    end
    field_out = field(i-1) + (field(i) - field(i-1)) .* (z_in - field_z(i-1))./(field_z(i) - field_z(i-1));
end

end