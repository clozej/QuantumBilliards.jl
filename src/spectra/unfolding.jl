

corner_correction(corner_angles) =  sum([(pi^2 - c^2)/(24*pi*c) for c in corner_angles])

weyl_law(A,L,k) = (A .* k.^2 .- L .* k)./(4*pi)
weyl_law(A,L,k,corner_angles) = weyl_law(A,L,k) .- corner_correction(corner_angles)