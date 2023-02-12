legmovement(k,a,b,c,offset) = a*cos(k*b*0.01*2*pi+offset)+c
Kp = [100;80;60]
Kd = [5;4;3]

function controller!(mechanism, k)
    angle21 = legmovement(k,paramcontainer[1][2],paramcontainer[1][1],paramcontainer[1][3],0)
    angle22 = legmovement(k,paramcontainer[1][2],paramcontainer[1][1],paramcontainer[1][3],pi)
    angle31 = legmovement(k,paramcontainer[1][4],paramcontainer[1][1],paramcontainer[1][5],-pi/2)
    angle32 = legmovement(k,paramcontainer[1][4],paramcontainer[1][1],paramcontainer[1][5],pi/2)

    for i=1:4
        θ1 = minimalCoordinates(mechanism, legjoints[(i-1)*3+1])[1]
        θ2 = minimalCoordinates(mechanism, legjoints[(i-1)*3+2])[1]
        θ3 = minimalCoordinates(mechanism, legjoints[(i-1)*3+3])[1]
        dθ1 = minimalVelocities(mechanism, legjoints[(i-1)*3+1])[1]
        dθ2 = minimalVelocities(mechanism, legjoints[(i-1)*3+2])[1]
        dθ3 = minimalVelocities(mechanism, legjoints[(i-1)*3+3])[1]

        if i == 1 || i == 4
            u1 = Kp[1]*(0-θ1) + Kd[1]*(0-dθ1)
            u2 = Kp[2]*(angle21-θ2) + Kd[2]*(0-dθ2)
            u3 = Kp[3]*(angle31-θ3) + Kd[3]*(0-dθ3)
        else
            u1 = Kp[1]*(0-θ1) + Kd[1]*(0-dθ1)
            u2 = Kp[2]*(angle22-θ2) + Kd[2]*(0-dθ2)
            u3 = Kp[3]*(angle32-θ3) + Kd[3]*(0-dθ3)
        end

        setForce!(mechanism, legjoints[(i-1)*3+1], SA[u1])
        setForce!(mechanism, legjoints[(i-1)*3+2], SA[u2])
        setForce!(mechanism, legjoints[(i-1)*3+3], SA[u3])
    end

    if trunk_body.state.xc[3]<0
        display("  upsidedown")
        setForce!(mechanism, legjoints[1], SA[10000.0]) # make simulation fail
    end
end