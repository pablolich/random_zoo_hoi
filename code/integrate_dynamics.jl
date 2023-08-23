using SciMLBase
using LinearAlgebra


function glv!(dy, y, s, W, t)
    Dy = Diagonal(y)
    dy =  Dy * (s - W*y)    
end