function shapeobject(cylinder::Cylinder)
    r,h = Tuple(cylinder.rh)
    GeometryTypes.Cylinder(Point(0.0,0.0,-h/2),Point(0.0,0.0,h/2), r)
end

function shapeobject(box::Box)
    x,y,z = Tuple(box.xyz)
    GeometryTypes.HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function shapeobject(box::CustomShape)
    x,y,z = Tuple(box.xyz)
    GeometryTypes.HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function visualize!(mechanism::Mechanism)
    vis = Visualizer()
    open(vis, Blink.Window())

    for shape in mechanism.shapes
        for id in shape.bodyids
            if id>=0
                # body = getbody(mechanism, id)
                visshape = shapeobject(shape)
                setobject!(vis["bundle/visshape"*string(id)], visshape, MeshPhongMaterial(color=shape.color))
            end
        end
    end

    framerate = Int64(round(1/mechanism.Δt))
    anim = MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate)

    for k=mechanism.steps
        MeshCat.atframe(anim, k) do
            for (id,body) in pairs(mechanism.bodies)
                settransform!(vis["bundle/visshape"*string(id)], compose(Translation(mechanism.storage.x[id][k]...),LinearMap(Quat(mechanism.storage.q[id][k]...))))
            end
        end
    end

    MeshCat.setanimation!(vis, anim)
    return
end


function convert_meshcat_to_video(;filename="video",input_path="util\\",output_path="util\\")
    # Saving MeshCat sequence as a video.
    meshcat_sequence_dir = joinpath(@__DIR__, "..", input_path)
    if filename==nothing
        filenames = readdir(meshcat_sequence_dir)
    else
        filenames = [filename * ".tar"]
    end
    for filename in filenames
        println("Converting " * filename * " to video." )
        video_dir = joinpath(@__DIR__, "..", output_path, filename[1:end-4] * ".mp4",)
        MeshCat.convert_frames_to_video(
            meshcat_sequence_dir * filename,
            video_dir,
            overwrite=true)
    end
    return
end
