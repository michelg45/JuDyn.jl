"""
    levicivita

        function computing the Levi-Civita  permutation symbol epsilon(i,j,k)

            Calling sequence:

                eps_ijk = levicivita(i,j,k)
"""
function levicivita(i::Int,j::Int,k::Int)

    if ( i > 3 || j >3 || k > 3 || i <= 0 || j <= 0 || k <= 0 )
        error("levicivita - incorrect values of i, j, k")
    end

    val = -1.0
    if ( i == j || j == k || k == i )
         val = 0.0
    elseif  (i < j && j < k) || (i < j && k < i) || (i > k && k > j)
        val = +1.0
    end
    return val
end