function parse_scalar(xel, name::String, T; default::String)
    if xel == nothing
        scalar = parse(T, default)
    else
        scalar = parse(T, attribute(xel, name))
    end

    return scalar
end

function parse_vector(xel, name::String, T; default::String)
    if xel == nothing || attribute(xel, name) == nothing
        vec = [parse(T, str) for str in split(default)]
    else
        vec = [parse(T, str) for str in split(attribute(xel, name))]
    end

    return vec
end

function parse_inertia_matrix(xinertia, T)
    if xinertia == nothing
        J = zeros(T, 3, 3)
    else
        ixx = parse_scalar(xinertia, "ixx", T, default = "0")
        ixy = parse_scalar(xinertia, "ixy", T, default = "0")
        ixz = parse_scalar(xinertia, "ixz", T, default = "0")
        iyy = parse_scalar(xinertia, "iyy", T, default = "0")
        iyz = parse_scalar(xinertia, "iyz", T, default = "0")
        izz = parse_scalar(xinertia, "izz", T, default = "0")
        J = [ixx ixy ixz; ixy iyy iyz; ixz iyz izz]
    end

    return J
end

function parse_pose(xpose, T)
    if xpose == nothing
        x, q = zeros(T, 3), Quaternion{T}()
    else
        x = parse_vector(xpose, "xyz", T, default = "0 0 0")
        rpy = parse_vector(xpose, "rpy", T, default = "0 0 0")
        q = Quaternion(RotZYX(rpy[3], rpy[2], rpy[1]))
    end

    return x, q
end

function parse_inertia(xinertial, T)
    if xinertial == nothing
        x = zeros(T, 3)
        q = Quaternion{T}()
        m = zero(T)
        J = zeros(T, 3, 3)
    else
        x, q = parse_pose(find_element(xinertial, "origin"), T)
        J = parse_inertia_matrix(find_element(xinertial, "inertia"), T)
        m = parse_scalar(find_element(xinertial, "mass"), "value", T, default = "0")
    end

    return x, q, m, J
end

function parse_xmaterial(xmaterial, T)
    if xmaterial == nothing
        color = RGBA(0.75, 0.75, 0.75)
    else
        colornode = find_element(xmaterial, "color")
        colorvec = parse_vector(colornode, "rgba", T, default = "0.5 0.5 0.5 1")
        color = RGBA(colorvec[1:3]...)
    end

    return color
end

function parse_shape(xvisual, T)
    if xvisual == nothing
        shape = nothing
    else
        xgeometry = find_element(xvisual, "geometry")
        @assert xgeometry != nothing

        color = parse_xmaterial(find_element(xvisual, "material"), T)
        x, q = parse_pose(find_element(xvisual, "origin"), T)

        shapenodes = LightXML.XMLElement[]
        for node in child_nodes(xgeometry)  # c is an instance of XMLNode
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

# TODO clean up relative shape to body frame
function parse_link(xlink, T)
    xvisual = find_element(xlink, "visual")

    x, q, m, J = parse_inertia(find_element(xlink, "inertial"), T)
    shape = parse_shape(find_element(xlink, "visual"), T)
    name = attribute(xlink, "name")

    if shape == nothing
        link = Body(m, J, name=name)
    else
        shape.m = m
        shape.J = J
        link = Body(shape, name=name)
    end

    link.x[1] = x 
    link.q[1] = q
    

    return link, shape
end

function parse_links(xlinks, T)
    ldict = Dict{String,AbstractBody{T}}()
    shapes = Shape{T}[]

    for xlink in xlinks
        link, shape = parse_link(xlink, T)
        ldict[link.name] = link
        shape != nothing && push!(shapes, shape)
    end

    return ldict, shapes
end

# TODO offset correct (e.g. for Planar)?
function parse_joint(xjoint, origin, plink, clink, T)
    joint_type = attribute(xjoint, "type")
    x, q = parse_pose(find_element(xjoint, "origin"), T)
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T, default = "1 0 0")
    p1 = x
    p2 = zeros(T, 3)
    name = attribute(xjoint, "name")
    
    # TODO limits for revolute joint?
    if joint_type == "revolute" || joint_type == "continuous"
        joint = EqualityConstraint(Revolute(plink, clink, p1, p2, axis, offset = q), name=name)
    elseif joint_type == "prismatic"
        joint = EqualityConstraint(Prismatic(plink, clink, p1, p2, axis, offset = q), name=name)
    elseif joint_type == "planar"
        joint = EqualityConstraint(Planar(plink, clink, p1, p2, axis), name=name)
    elseif joint_type == "fixed"
        joint = EqualityConstraint(Fixed(plink, clink, p1, p2, offset = q), name=name)
    elseif joint_type == "floating" # Floating relative to non-origin link not supported
        @assert origin == plink
        joint = EqualityConstraint(OriginConnection(origin, clink), name=name)
    end

    return joint
end

function parse_joints(xjoints, ldict, T, floating)
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

function parse_urdf(filename, ::Type{T}, floating) where T
    xdoc = LightXML.parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xlinks = get_elements_by_tagname(xroot, "link")
    xjoints = get_elements_by_tagname(xroot, "joint")

    ldict, shapes = parse_links(xlinks, T)
    origin, links, joints = parse_joints(xjoints, ldict, T, floating)

    free(xdoc)

    return origin, links, joints, shapes
end