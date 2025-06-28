"""
    Lmesh(maxSQlen::Float64)

# Arguments
 - maxSQlen::Float64 : The maximum side size of the square with which you want to create the L mesh. (This is the max length, if it doesn't 
                        divide the domain equally the function will take nearest value below maxSQlen which creates an even mesh)

# Returns
- go : Global Ordering (go) is a matrix with 3 other martices (the subdomains) inside it, those 3 matrices contain the mesh of each subdomain 
"""
function Lmesh(maxSQlen::Float64)

    domain = [-1 0 1;] # Domain

    # Func input
    domStep = abs(domain[2] - domain[1]);
    maxSQlen = min(maxSQlen, domStep); # Maximum square length (not greater than 1 so that we dont overwrite the domain)
    maxSQ = div(domStep, maxSQlen, RoundUp); # Maximum small squares inside the big square based on maxSQlen
    step1 =  domStep/maxSQ; # Final closest step so that all small squares have equal width <= maxSQlen
    
    #Domain properties
    n = length(domain);
    domMid = (n+1)/2;
    
    # ---- Square 1 & 2 Mesh --------------------------------
    
    x1 = range(domain[1], domain[Integer(domMid)], step = step1) |> collect; # Nodes vector in row 1
    rowNodes = length(x1);
    rowElem = length(x1) - 1; # Elements in 1 row
    
    # Global ordering of Square 1 & 2
    
    # Faster way than for loop, but same complexity, O(n^2):
    go1 = (rowNodes-1:-1:0).*(2*rowNodes-1).+(1:rowNodes)'; # Only use matrix operations to create the mesh
    lastNode = go1[end,end]; # Temp to hold the last node so that we can construct the next matrix
    
    go2 = (rowNodes-1:-1:0).*(2*rowNodes-1).+(lastNode:lastNode+rowNodes-1)';
    lastNode = go2[1,end];
    
    # ---- Square 3 Mesh --------------------------------
    
    go3 = zeros(Int, rowElem+1, rowElem+1);
    go3[end,:] = go1[1,:];
    go3[1:end-1, :] = (rowNodes-2:-1:0).*(rowNodes).+(lastNode+1:lastNode+rowNodes)';
    
    go = [go1, go2, go3]; # Use go[1] to get go1, go[2] to get go2 etc.

    return go
end