function shapeobject(cylinder::Cylinder)
    r,h = Tuple(cylinder.rh)
    GeometryTypes.Cylinder(Point(0.0,0.0,-h/2),Point(0.0,0.0,h/2), r)
end

function shapeobject(box::Box)
    x,y,z = Tuple(box.xyz)
    GeometryTypes.HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function visualize(mechanism::Mechanism,shapes)
    vis = Visualizer()
    open(vis, Blink.Window())
    # open(vis)

    for body in mechanism.bodies
        for shape in shapes
            for shapebody in shape.bodies
                if shapebody == body
                    visbody = shapeobject(shape)
                    setobject!(vis["bundle/visbody"*string(body.id)], visbody, MeshPhongMaterial(color=shape.color))
                    break
                end
            end
        end
    end

    framerate = Int64(round(1/mechanism.dt))
    anim = MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate)

    for k=mechanism.steps
        MeshCat.atframe(anim, vis, k) do frame
            for (id,body) in pairs(mechanism.bodies)
                settransform!(vis["bundle/visbody"*string(id)], compose(Translation(mechanism.storage.x[id][k]...),LinearMap(Quat(mechanism.storage.q[id][k]...))))
            end
        end
    end

    MeshCat.setanimation!(vis, anim);
end
