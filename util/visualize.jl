function GeometryTypes.Cylinder(box::Box)
    Cylinder(Point(0.0,0.0,-box.lwh[3]/2),Point(0.0,0.0,box.lwh[3]/2), box.lwh[1]/2)
end

function visualize(robot::Robot,shapes)
    vis = Visualizer()
    open(vis, Blink.Window())

    vislinks = Vector{GeometryPrimitive}(undef,0)
    for link in robot.links
        for shape in shapes
            for id in shape.linkids
                if id == link.id
                    vislink = Cylinder(shape)
                    setobject!(vis["bundle/vislink"*string(id)], vislink, MeshPhongMaterial(color=shape.color))
                    push!(vislinks,Cylinder(shape))
                    break
                end
            end
        end
    end

    framerate = Int64(round(1/robot.dt))
    anim = MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate)

    for k=robot.steps
        # Set the pose of the car bundle.
        # MeshCat.atframe(anim, vis, k) do frame
        #     settransform!(vis["bundle/vislink1"], compose(Translation(bot.storage.x[1][k]...),LinearMap(Quat(bot.storage.q[1][k]...))))
        #     settransform!(vis["bundle/vislink2"], compose(Translation(bot.storage.x[2][k]...),LinearMap(Quat(bot.storage.q[2][k]...))))
        #     settransform!(vis["bundle/vislink3"], compose(Translation(bot.storage.x[3][k]...),LinearMap(Quat(bot.storage.q[3][k]...))))
        #     settransform!(vis["bundle/vislink4"], compose(Translation(bot.storage.x[4][k]...),LinearMap(Quat(bot.storage.q[4][k]...))))
        # end
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
