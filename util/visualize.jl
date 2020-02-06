function shapeobject(cylinder::Cylinder)
    r,h = Tuple(cylinder.rh)
    GeometryTypes.Cylinder(Point(0.0,0.0,-h/2),Point(0.0,0.0,h/2), r)
end

function shapeobject(box::Box)
    x,y,z = Tuple(box.xyz)
    GeometryTypes.HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function visualize(robot::Robot,shapes)
    vis = Visualizer()
    open(vis, Blink.Window())
    # open(vis)

    for link in robot.links
        for shape in shapes
            for shapelink in shape.links
                if shapelink == link
                    vislink = shapeobject(shape)
                    setobject!(vis["bundle/vislink"*string(link.id)], vislink, MeshPhongMaterial(color=shape.color))
                    break
                end
            end
        end
    end

    framerate = Int64(round(1/robot.dt))
    anim = MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate)

    for k=robot.steps
        MeshCat.atframe(anim, vis, k) do frame
            for (id,link) in pairs(robot.links)
                settransform!(vis["bundle/vislink"*string(id)], compose(Translation(robot.storage.x[id][k]...),LinearMap(Quat(robot.storage.q[id][k]...))))
            end
        end
    end

    MeshCat.setanimation!(vis, anim);
end
