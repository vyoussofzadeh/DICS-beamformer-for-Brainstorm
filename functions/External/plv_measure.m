function plv = plv_measure(data)
[ nc, ns, nt ] = size( data );
ndat = data ./ abs( data );
plv = zeros( nc, nc, nt );
for t = 1: nt
    plv( :, :, t ) = abs( ndat( :, :, t ) * ndat( :, :, t )' ) / ns;
end
end