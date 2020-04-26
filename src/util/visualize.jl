function shapeobject(box::Box)
    x,y,z = Tuple(box.xyz)
    return GeometryTypes.HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function shapeobject(cylinder::Cylinder)
    r,h = Tuple(cylinder.rh)
    return GeometryTypes.Cylinder(Point(0.0,0.0,-h/2),Point(0.0,0.0,h/2), r)
end

function shapeobject(sphere::Sphere)
    r = sphere.r
    return GeometryTypes.Sphere(Point(0.0,0.0,0.0), r)
end

function shapeobject(mesh::Mesh)
    return shape = load(mesh.path, GLUVMesh)
end

function visualize!(mechanism::Mechanism)
    vis = Visualizer()
    open(vis, Blink.Window())

    storage = deepcopy(mechanism.storage)
    oid = mechanism.origin.id
    oshapeind = 0

    for (ind,shape) in enumerate(mechanism.shapes)
        for id in shape.bodyids
            if id >= 0
                for i in mechanism.steps
                    storage.x[id][i] += vrotate(shape.xoff, storage.q[id][i])
                    storage.q[id][i] *= shape.qoff
                end
                visshape = shapeobject(shape)
                setobject!(vis["bundle/visshape"*string(id)], visshape, MeshPhongMaterial(color=shape.color))
            else
                @assert id == oid
                oshapeind = ind
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
                settransform!(vis["bundle/visshape"*string(id)], compose(Translation((storage.x[id][k])...),LinearMap(Quat((storage.q[id][k])...))))
            end
            if oshapeind > 0
                shape = mechanism.shapes[oshapeind]
                settransform!(vis["bundle/visshape"*string(oid)], compose(Translation((shape.xoff)...),LinearMap(Quat((shape.qoff)...))))
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
