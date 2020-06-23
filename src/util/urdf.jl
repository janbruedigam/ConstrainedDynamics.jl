unsafeattribute(::Nothing, ::Core.AbstractString) = nothing
unsafeattribute(x::LightXML.XMLElement, name::Core.AbstractString) = attribute(x, name)

function parse_scalar(xel, name::String, T; default::Union{String,Nothing}=nothing)
    scalarstr = unsafeattribute(xel, name)
    if scalarstr === nothing 
        if default === nothing
            @error "no parsable scalar found"
        else
            scalarstr = default
        end
    end

    return parse(T, scalarstr)
end

function parse_vector(xel, name::String, T; default::Union{String,Nothing}=nothing)
    vectorstr = unsafeattribute(xel, name)
    if vectorstr === nothing 
        if default === nothing
            @error "no parsable vector found"
        else
            vectorstr = default
        end
    end

    return parse.(T,split(vectorstr))
end

function parse_inertiamatrix(xinertia, T)
    if xinertia === nothing
        J = zeros(T, 3, 3)
    else
        ixx = parse_scalar(xinertia, "ixx", T)
        ixy = parse_scalar(xinertia, "ixy", T, default = "0")
        ixz = parse_scalar(xinertia, "ixz", T, default = "0")
        iyy = parse_scalar(xinertia, "iyy", T)
        iyz = parse_scalar(xinertia, "iyz", T, default = "0")
        izz = parse_scalar(xinertia, "izz", T)
        J = [ixx ixy ixz; ixy iyy iyz; ixz iyz izz]
    end

    return J
end

function parse_pose(xpose, T)
    if xpose === nothing
        x, q = zeros(T, 3), one(UnitQuaternion{T})
    else
        x = parse_vector(xpose, "xyz", T, default = "0 0 0")
        rpy = parse_vector(xpose, "rpy", T, default = "0 0 0")
        q = UnitQuaternion(RotZYX(rpy[3], rpy[2], rpy[1]))
    end

    return x, q
end

function parse_inertia(xinertial, T)
    if xinertial === nothing
        x = zeros(T, 3)
        q = one(UnitQuaternion{T})
        m = zero(T)
        J = zeros(T, 3, 3)
    else
        x, q = parse_pose(find_element(xinertial, "origin"), T)
        J = parse_inertiamatrix(find_element(xinertial, "inertia"), T)
        m = parse_scalar(find_element(xinertial, "mass"), "value", T, default = "0")
    end

    return x, q, m, J
end

function parse_robotmaterials(xroot, T)
    xmaterials = get_elements_by_tagname(xroot, "material")

    mdict = Dict{String,Vector{T}}()

    for xmaterial in xmaterials
        name = attribute(xmaterial, "name")
        colornode = find_element(xmaterial, "color")
        colorvec = parse_vector(colornode, "rgba", T)
        mdict[name] = colorvec
    end

    return mdict
end

function parse_xmaterial(xmaterial, materialdict, T)
    if xmaterial === nothing
        color = RGBA(0.75, 0.75, 0.75)
    else
        name = unsafeattribute(xmaterial, "name")
        if haskey(materialdict, name)
            colorvec = materialdict[name]
        else
            colorvec = [0.75; 0.75; 0.75; 1.0]
        end
        colornode = find_element(xmaterial, "color")
        colorvec = parse_vector(colornode, "rgba", T, default = string(colorvec[1]," ",colorvec[2]," ",colorvec[3]," ",colorvec[4]))

        color = RGBA(colorvec...)
    end

    return color
end

function parse_shape(xvisual, materialdict, T)
    if xvisual === nothing
        shape = nothing
    else
        xgeometry = find_element(xvisual, "geometry")
        @assert xgeometry !== nothing

        color = parse_xmaterial(find_element(xvisual, "material"), materialdict, T)
        x, q = parse_pose(find_element(xvisual, "origin"), T)

        shapenodes = LightXML.XMLElement[]
        for node in child_nodes(xgeometry)  # node is an instance of XMLNode
            if is_elementnode(node)
                push!(shapenodes, XMLElement(node))
            end
        end

        if length(shapenodes) == 0
            shape = nothing
        else 
            if length(shapenodes) > 1
                @info "Multiple geometries."
            end

            shapenode = shapenodes[1]
            if name(shapenode) == "box"
                xyz = parse_vector(shapenode, "size", T, default = "1 1 1")
                shape = Box(xyz..., zero(T), color = color, xoff = x, qoff = q)
            elseif name(shapenode) == "cylinder"
                r = parse_scalar(shapenode, "radius", T, default = "0.5")
                l = parse_scalar(shapenode, "length", T, default = "1")
                shape = Cylinder(r, l, zero(T), color = color, xoff = x, qoff = q)
            elseif name(shapenode) == "sphere"
                r = parse_scalar(shapenode, "radius", T, default = "0.5")
                shape = Sphere(r, zero(T), color = color, xoff = x, qoff = q)
            elseif name(shapenode) == "mesh"
                path = attribute(shapenode, "filename")
                shape = Mesh(path, zero(T), zeros(T, 3, 3), color = color, xoff = x, qoff = q)
            else
                @info "Unkown geometry."
                shape = nothing
            end
        end

        
    end

    return shape
end

function parse_link(xlink, materialdict, T)
    xvisual = find_element(xlink, "visual")

    x, q, m, J = parse_inertia(find_element(xlink, "inertial"), T)
    shape = parse_shape(find_element(xlink, "visual"), materialdict, T)
    name = attribute(xlink, "name")

    if shape === nothing
        link = Body(m, J, name=name)
    else
        shape.m = m
        shape.J = J
        link = Body(shape, name=name)
    end

    link.state.xc = x 
    link.state.qc = q
    

    return link, shape
end

function parse_links(xlinks, materialdict, T)
    ldict = Dict{String,AbstractBody{T}}()
    shapes = Shape{T}[]

    for xlink in xlinks
        link, shape = parse_link(xlink, materialdict, T)
        ldict[link.name] = link
        shape !== nothing && push!(shapes, shape)
    end

    return ldict, shapes
end

# TODO offset correct (e.g. for Planar)?
function parse_joint(xjoint, origin, plink, clink, T)
    jointtype = attribute(xjoint, "type")
    x, q = parse_pose(find_element(xjoint, "origin"), T)
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T, default = "1 0 0")
    p1 = x
    p2 = zeros(T, 3)
    name = attribute(xjoint, "name")
    
    # TODO limits for revolute joint?
    if jointtype == "revolute" || jointtype == "continuous"
        joint = EqualityConstraint(Revolute(plink, clink, axis; p1=p1, p2=p2, qoffset = q), name=name)
    elseif jointtype == "prismatic"
        joint = EqualityConstraint(Prismatic(plink, clink, axis; p1=p1, p2=p2, qoffset = q), name=name)
    elseif jointtype == "planar"
        joint = EqualityConstraint(Planar(plink, clink, axis; p1=p1, p2=p2), name=name)
    elseif jointtype == "fixed"
        joint = EqualityConstraint(Fixed(plink, clink; p1=p1, p2=p2, qoffset = q), name=name)
    elseif jointtype == "floating" # Floating relative to non-origin link not supported
        @assert origin == plink
        joint = EqualityConstraint(OriginConnection(origin, clink), name=name)
    end

    return joint
end

function parse_joints(xjoints, ldict, floating, T)
    origins = Origin{T}[]
    links = Body{T}[]
    joints = EqualityConstraint{T}[]
    floatingname = ""

    for name in keys(ldict)
        child = false
        for xjoint in xjoints
            xchild = find_element(xjoint, "child")
            childname = attribute(xchild, "link")
            if childname == name
                child = true
                break
            end
        end
        if child
            push!(links, ldict[name])
        else
            origin = Origin{T}()
            if floating # keep current link and create new origin
                push!(links, ldict[name])
                floatingname = name
            else # make current link origin
                origin.id = ldict[name].id
                origin.name = name
            end
            push!(origins, origin)
        end
    end

    @assert length(origins) == 1 "Multiple origins"
    origin = origins[1]

    for xjoint in xjoints
        xplink = find_element(xjoint, "parent")
        plink = ldict[attribute(xplink, "link")]
        if plink.id == origin.id
            plink = origin
        end
        xclink = find_element(xjoint, "child")
        clink = ldict[attribute(xclink, "link")]

        joint = parse_joint(xjoint, origin, plink, clink, T)
        push!(joints, joint)
    end

    if floating
        originjoint = EqualityConstraint(OriginConnection(origin, ldict[floatingname]),name="autoorigincon")
        push!(joints, originjoint)
    end

    return origin, links, joints
end

function parse_urdf(filename, floating, ::Type{T}) where T
    xdoc = LightXML.parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xlinks = get_elements_by_tagname(xroot, "link")
    xjoints = get_elements_by_tagname(xroot, "joint")

    materialdict = parse_robotmaterials(xroot, T)
    ldict, shapes = parse_links(xlinks, materialdict, T)
    origin, links, joints = parse_joints(xjoints, ldict, floating, T)

    free(xdoc)

    return origin, links, joints, shapes
end