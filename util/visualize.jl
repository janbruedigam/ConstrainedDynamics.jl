# function GeometryTypes.Cylinder(box::Box)
#     Cylinder(Point(0.0,0.0,-box.lwh[3]/2),Point(0.0,0.0,box.lwh[3]/2), box.lwh[1]/2)
# end

function GeometryTypes.HyperRectangle(box::Box)
    x,y,z = Tuple(box.xyz)
    HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function visualize(robot::Robot,shapes)
    vis = Visualizer()
    open(vis, Blink.Window())

    for link in robot.links
        for shape in shapes
            for id in shape.linkids
                if id == link.id
                    vislink = HyperRectangle(shape)
                    setobject!(vis["bundle/vislink"*string(id)], vislink, MeshPhongMaterial(color=shape.color))
                    break
                end
            end
        end
    end

    framerate = Int64(round(1/robot.dt))
    anim = MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate)

    for k=robot.steps
        MeshCat.atframe(anim, vis, k) do frame
            for link in robot.links
                id = link.id
                ind = robot.ldict[id]
                settransform!(vis["bundle/vislink"*string(id)], compose(Translation(robot.storage.x[ind][k]...),LinearMap(Quat(robot.storage.q[ind][k]...))))
            end
        end
    end

    MeshCat.setanimation!(vis, anim);
end
