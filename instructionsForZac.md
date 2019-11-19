Relevant files:
test.jl: Run this file for the full simulation
simulation.jl: Enter newton's method from this file
newton.jl: Newton's method, line search, and iterative refinement

test.jl:
Lines 55 to 83
Select the mechanism you want to simulate by commenting in/out the respective lines
    - eg uncomment lines 57 and 58 for the closed chain 4 bar mechanism
    - or uncomment lines 61 and 62 for the open chain 4 bar mechanism

simulation.jl:
Line 38
Enter newton's method here

To debug, comment out lines 78 and 79,
    run         include("test.jl") once
    then run    Juno.@enter sim!(simul,debug=false,disp=true)
    step into   sol = newton(robot,warning=debug)[1]      (line 38 in simulation.jl)
