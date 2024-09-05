#include("../abstracttypes.jl")

"""
Check if two numbers are equal within specified tolerances.

# Arguments
- `x`: The first number.
- `dx`: The tolerance for the first number.
- `y`: The second number.
- `dy`: The tolerance for the second number.

# Returns
- `Bool`: `true` if the intervals `[x - dx, x + dx]` and `[y - dy, y + dy]` overlap, `false` otherwise.

# Logic
The function checks if the intervals `[x - dx, x + dx]` and `[y - dy, y + dy]` overlap.
"""
function is_equal(x::T, dx::T, y::T, dy::T) :: Bool where {T<:Real}
    # Define the intervals
    x_lower = x - dx
    x_upper = x + dx
    y_lower = y - dy
    y_upper = y + dy
    # Check if the intervals overlap
    return max(x_lower, y_lower) <= min(x_upper, y_upper)
end

"""
Match and merge two sorted lists of wavenumbers with associated tensions.

# Arguments
- `ks_l`: A sorted vector of wavenumbers from the left set.
- `ts_l`: A vector of associated tensions corresponding to `ks_l`.
- `ks_r`: A sorted vector of wavenumbers from the right set.
- `ts_r`: A vector of associated tensions corresponding to `ks_r`.

# Returns
- `ks`: A vector of merged wavenumbers, taking into account overlaps.
- `ts`: A vector of corresponding tensions for the merged wavenumbers.
- `control`: A vector of booleans indicating whether a wavenumber was matched (`true`) or taken from a single list (`false`).

# Logic
- The function iterates through both sorted lists of wavenumbers, `ks_l` and `ks_r`, and checks for matches using the `is_equal` function.
- If a match is found (i.e., the intervals overlap), the wavenumber with the smaller tension is kept, and both indices are incremented. To the control vector we also add `true` to indicate that an overlap for that k was found.
- If no match is found, the function adds the smaller wavenumber to the result, and the index of that list is incremented.To the control vector we also add `false` to indicate that no overlap for that k was found.
- The process continues until all wavenumbers are processed, ensuring that the resulting vectors `ks` and `ts` are sorted and correctly merged.

"""
function match_wavenumbers(ks_l::Vector{T},ts_l::Vector{T},ks_r::Vector{T},ts_r::Vector{T}) where {T<:Real}
    #vectors ks_l and_ks_r must be sorted. Sanity check?
    p_l = sortperm(ks_l) # get the sorted indexes for left
    p_r = sortperm(ks_r) # get the sorted indexes for right
    # Order them by indexes of ks left and right respectively
    ks_l = ks_l[p_l] 
    ts_l = ts_l[p_l]
    ks_r = ks_r[p_r]
    ts_r = ts_r[p_r]
    
    #vectors ks_l and_ks_r must be sorted
    i = j = 1 #counting index
    control = Vector{Bool}()#control bits
    ks = Vector{eltype(ks_l)}()#final wavenumbers
    ts = Vector{eltype(ts_l)}()#final tensions
    while i <= length(ks_l) && j <= length(ks_r)
        x, dx = ks_l[i], ts_l[i]
        y, dy = ks_r[j], ts_r[j]
        if  is_equal(x,dx,y,dy) #check equality with errorbars
            i += 1 
            j += 1
            if dx < dy
                push!(ks, x)
                push!(ts, dx)
                push!(control, true)
            else
                push!(ks, y)
                push!(ts, dy)
                push!(control, true)
            end
        elseif x < y
            i += 1
            push!(ks, x)
            push!(ts, dx)
            push!(control, false)
        else 
            j += 1
            push!(ks, y)
            push!(ts, dy)
            push!(control, false)
        end
    end
    return ks, ts, control 
end

# Skupaj združi majhne spektre, za vsak dk zmerga skupaj spekter. Pogledaš overlaping intervale in vsak level doibš vakrat, enkrat z leve in desne. Za vsak level moraš pogledati če je v level ali desnem hkrati. Imaš ja ali ne kriteriij ki vzame bolj natančnega

# Dinamični kriteriij za določitev dk z s pomočjo Weylove formule. Damo nekaj srednjih razmikov notri in kako narašča. Barbara Dietz ima v članku kriterij kako izbirati. To bo praktičen kriterij

# Določi najbolj natančne nivoje. Tensions morajo biti 1e-8, včasih sprejmeš tudi 1e-4
# Napaka mora biti precej manjša kot srednji razmik med nivoji -> končen precision in kako se to zgodi
# Velikost baze imaš failsafe na koncu zaradi velikosti baze. Benchmark narediš da vidiš velikost matrike pred in po diagonalizaciji -> lahko zmanjšaš skalirni faktor če je nuja

# Imamo cel spekter in želim preverit k0 in pogledaš v starem spektru pogledaš in zmergaš noter -> DataFrames

# !!! Težava je IntervalArithmetic in to je treba odpraviti (druga pa k0 problem preblizu in ga odreže). k-ji so do tensiona natančni! -> težava je v is_equal!
# Za mushroom parametrizacija je Gauss-Legendre nodes

# Radially integrated and Angularly integrated -> tukaj bi lahko še članke dodali
# Kako vidiš da imaš napako -> še enkrat naokrog! -> mora biti relativo preprosta rešitev (upoštevaj da imajo napake in pazi na interval prekrivanja)

# Unit tests -> kakšnen nivoje moraš odbit, vse funkcije kličeš da vidiš ali dobiš vse
# Vedno pazi na tope kote -> sploh pri trikotniku

# !!! Shift angle -pi/2 pri bazi za half mushroom (tukaj je to izredno pomembno, zadnjo kot je 3pi/2) -> nujno preveri 
# Gap ratios -> naredi distribucije za to (debeli repi moramo nomralizacijio pazit!) -> potenčen rep in imapš končen spekter in boš samo od nekega gapa imel podatke in če greš dovolj nizko bo normalizacija nenatannčna. Potrebno je paziti in izračunati dodaten integral za rep. Koliko verjetnosti v repu moraš odštet od repa da dobiš normalizacijo

"""
Merge and resolve overlapping wavenumbers from two sets, modifying the first set in place.

# Arguments
- `k_left`: The vector of wavenumbers from the left set, which will be modified in place.
- `ten_left`: The vector of tensions associated with `k_left`, modified in place.
- `k_right`: The vector of wavenumbers from the right set.
- `ten_right`: The vector of tensions associated with `k_right`.
- `control_left`: A control vector corresponding to `k_left`, indicating if wavenumbers were matched, modified in place.
- `kl`: The lower bound of the overlap interval.
- `kr`: The upper bound of the overlap interval.
- `tol`: Tolerance for defining the overlap interval, default is `1e-3`.

# Returns
- Modifies `k_left`, `ten_left`, and `control_left` in place to reflect the merged wavenumbers and tensions.

# Logic
1. **Identify Overlapping Intervals**: 
   - Determine which elements in `k_left` and `k_right` lie within the overlap interval `[kl - tol, kr + tol]`.
   
2. **Match Wavenumbers**:
   - Use `match_wavenumbers` to identify matching wavenumbers between the overlapping segments of `k_left` and `k_right`, resolving them into a unified list of wavenumbers and tensions.

3. **Modify `k_left` In-Place**:
   - Remove the overlapping wavenumbers and their associated tensions from `k_left`, `ten_left`, and `control_left`.
   - Add the matched wavenumbers and tensions into `k_left` and `ten_left`, updating `control_left` to indicate which wavenumbers were successfully matched.

4. **Append Remaining Elements**:
   - After handling the overlap, any wavenumbers and tensions from `k_right` that are outside the overlap region are added to the end of `k_left` and `ten_left`. These new additions are marked in `control_left` as unmatched (`false`), indicating they were not part of the original overlapping region.

# Reason for `control_left` Modification:
- The boolean values in `control_left` are updated to reflect whether each wavenumber was matched with a corresponding value from the other set. This helps track which wavenumbers were directly matched and which were added as unmatched elements.
"""
function overlap_and_merge!(k_left::Vector{T}, ten_left::Vector{T}, k_right::Vector{T}, ten_right::Vector{T}, control_left::Vector{Bool}, kl::T, kr::T; tol=1e-3) where {T<:Real}
    #find overlaps in interval [k1,k2]
    idx_l = k_left .> (kl-tol) .&& k_left .< (kr+tol)
    idx_r = k_right .> (kl-tol) .&& k_right .< (kr+tol)
    
    ks_l,ts_l,ks_r,ts_r = k_left[idx_l], ten_left[idx_l], k_right[idx_r], ten_right[idx_r]
    #check if wavnumbers match in overlap interval
    ks, ts, control = match_wavenumbers(ks_l,ts_l,ks_r,ts_r)
    #println("left: $ks_l")
    #println("right: $ks_r")
    #println("overlaping: $ks")
    #i_l = idx_l[1]
    #i_r = idx_r[end]+1
    deleteat!(k_left, idx_l)
    append!(k_left, ks)
    deleteat!(ten_left, idx_l)
    append!(ten_left, ts)
    deleteat!(control_left, idx_l)
    append!(control_left, control)

    idx_last = findlast(idx_r) + 1
    append!(k_left, k_right[idx_last:end])
    append!(ten_left, ten_right[idx_last:end])
    append!(control_left, [false for i in idx_last:length(k_right)])
end

"""
Compute the spectrum of wavenumbers within a specified range. This is a high level wrapper that should be used in palce of it's lower level constituents

# Arguments
- `solver::AbsSolver`: The solver instance used for solving the eigenproblem.
- `basis::AbsBasis`: The basis functions used in the problem.
- `pts::AbsPoints`: The evaluated points on the boundary or domain of the problem.
- `k1`: The lower bound of the wavenumber range.
- `k2`: The upper bound of the wavenumber range.
- `dk`: The step size for incrementing the wavenumber.
- `tol`: Tolerance for resolving overlaps between wavenumbers, default is `1e-4`.
- `plot_info`: A boolean flag indicating whether to plot or output additional information (currently not utilized), default is `false`.

# Returns
- `k_res`: A vector of computed wavenumbers across the range `[k1, k2]`.
- `ten_res`: A vector of corresponding tensions for each wavenumber.
- `control`: A vector of booleans indicating whether each wavenumber was matched during the merging process (left & right).

# Logic
1. **Initial Computation**: 
   - Start by solving the eigenproblem at the initial wavenumber `k1`, storing the resulting wavenumbers and tensions in `k_res` and `ten_res`. Initialize a control vector to track matched wavenumbers.

2. **Iterative Solving**:
   - Incrementally step through the wavenumber range from `k1` to `k2` by `dk`, solving the eigenproblem at each step. For each new set of wavenumbers and tensions, use `overlap_and_merge!` to merge them into the results, resolving overlaps and updating the control vector.

3. **Overlap Resolution**:
   - The `overlap_and_merge!` function ensures that wavenumbers and tensions from consecutive steps are appropriately merged, accounting for overlaps within the specified tolerance `tol`.

4. **Final Output**:
   - After iterating through the entire range, return the computed wavenumbers, their corresponding tensions, and the control vector.

# Reason for Control Vector:
- The control vector is used to track whether a wavenumber was present in both intervals (left & right) in the overlap_and_merge function.
"""
function compute_spectrum(solver::AbsSolver, basis::AbsBasis, pts::AbsPoints,k1,k2,dk;tol=1e-4, plot_info=false)
    k0 = k1
    #initial computation
    k_res, ten_res = solve(solver, basis, pts, k0, dk+tol)
    control = [false for i in 1:length(k_res)]
    while k0 < k2
        k0 += dk
        k_new, ten_new = solve(solver, basis, pts, k0, dk+tol)
        overlap_and_merge!(k_res, ten_res, k_new, ten_new, control, k0-dk, k0; tol=tol)

    end
    return k_res, ten_res, control
end