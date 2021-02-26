### Before parsing

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
                shape = Box(xyz..., zero(T), color = color, xoffset = x, qoffset = q)
            elseif name(shapenode) == "cylinder"
                r = parse_scalar(shapenode, "radius", T, default = "0.5")
                l = parse_scalar(shapenode, "length", T, default = "1")
                shape = Cylinder(r, l, zero(T), color = color, xoffset = x, qoffset = q)
            elseif name(shapenode) == "sphere"
                r = parse_scalar(shapenode, "radius", T, default = "0.5")
                shape = Sphere(r, zero(T), color = color, xoffset = x, qoffset = q)
            elseif name(shapenode) == "mesh"
                path = attribute(shapenode, "filename")
                scale = parse_vector(shapenode, "scale", T, default = "1 1 1")
                shape = Mesh(path, zero(T), zeros(T, 3, 3), scale=scale, color = color, xoffset = x, qoffset = q)
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
        link = shape
        link.m = m
        link.J = J
        link.name = name
    end

    link.state.xc = x 
    link.state.qc = q
    

    return link
end

function parse_links(xlinks, materialdict, T)
    ldict = Dict{String,AbstractBody{T}}()

    for xlink in xlinks
        link = parse_link(xlink, materialdict, T)
        ldict[link.name] = link
    end

    return ldict
end

function parse_joint(xjoint, plink, clink, T)
    jointtype = attribute(xjoint, "type")
    x, q = parse_pose(find_element(xjoint, "origin"), T)
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T, default = "1 0 0")
    p1 = x
    name = attribute(xjoint, "name")
    
    # TODO limits for revolute joint?
    if jointtype == "revolute" || jointtype == "continuous"
        joint = EqualityConstraint(Revolute(plink, clink, axis; p1=p1, qoffset = q), name=name)
    elseif jointtype == "prismatic"
        joint = EqualityConstraint(Prismatic(plink, clink, axis; p1=p1, qoffset = q), name=name)
    elseif jointtype == "planar"
        joint = EqualityConstraint(Planar(plink, clink, axis; p1=p1, qoffset = q), name=name)
    elseif jointtype == "fixed"
        joint = EqualityConstraint(Fixed(plink, clink; p1=p1, qoffset = q), name=name)
    elseif jointtype == "floating"
        joint = EqualityConstraint(Floating(plink, clink), name=name)
    else
        @error "Unknown joint type"
    end

    return joint
end

function parse_loop_joint(xjoint, link1, link2, T)
    find_element(xjoint, "link1")
    find_element(xjoint, "link2")

    jointtype = attribute(xjoint, "type")
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T, default = "1 0 0")
    x1, q1 = parse_pose(find_element(xjoint, "link1"), T)
    x2, q2 = parse_pose(find_element(xjoint, "link2"), T) # The orientation q2 of the second body is ignored because it is determined by the mechanism's structure
    p1 = x1
    p2 = x2
    name = attribute(xjoint, "name")
    
    # TODO limits for revolute joint?
    if jointtype == "revolute" || jointtype == "continuous"
        joint = EqualityConstraint(Revolute(link1, link2, axis; p1=p1, p2=p2, qoffset = q1), name=name)
    elseif jointtype == "prismatic"
        joint = EqualityConstraint(Prismatic(link1, link2, axis; p1=p1, p2=p2, qoffset = q1), name=name)
    elseif jointtype == "planar"
        joint = EqualityConstraint(Planar(link1, link2, axis; p1=p1, p2=p2, qoffset = q1), name=name)
    elseif jointtype == "fixed"
        joint = EqualityConstraint(Fixed(link1, link2; p1=p1, p2=p2, qoffset = q1), name=name)
    elseif jointtype == "floating"
        joint = EqualityConstraint(Floating(link1, link2), name=name)
    else
        @error "Unknown joint type"
    end

    return joint
end

function parse_joints(xjoints, ldict, floating, T)
    origins = Origin{T}[]
    links = Body{T}[]
    joints = EqualityConstraint{T}[]
    floatingname = ""

    for name in keys(ldict)
        childflag = false
        for xjoint in xjoints
            xchild = find_element(xjoint, "child")
            childname = attribute(xchild, "link")
            if childname == name
                childflag = true
                break
            end
        end
        if childflag
            push!(links, ldict[name])
        else
            origin = Origin{T}()
            if floating # keep current link and create new origin
                push!(links, ldict[name])
                floatingname = name
            else # make current link origin
                origin = Origin(ldict[name])
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

        joint = parse_joint(xjoint, plink, clink, T)
        push!(joints, joint)
    end

    if floating
        originjoint = EqualityConstraint(Floating(origin, ldict[floatingname]), name="auto_generated_floating_joint")
        push!(joints, originjoint)
    end

    return origin, links, joints
end

# TODO This might be missing the detection of a direct loop, i.e. only two links connected by two joints
# TODO Also only works for a single loop closure in a cycle (so no ladders)
function parse_loop_joints(xloopjoints, origin, joints, ldict, T)
    loopjoints = EqualityConstraint{T}[]


    for xloopjoint in xloopjoints
        xlink1 = find_element(xloopjoint, "link1")
        xlink2 = find_element(xloopjoint, "link2")
        link1 = ldict[attribute(xlink1, "link")]
        link2 = ldict[attribute(xlink2, "link")]

        predlist = Tuple{Int64,Int64}[]
        jointlist = [(joints[i].id,joints[i].parentid,joints[i].childids) for i=1:length(joints)]
        linkid = link1.id

        while true # create list of predecessor joints and parent links for link1
            for (i,jointdata) in enumerate(jointlist)
                if linkid ∈ jointdata[3]
                    push!(predlist,(jointdata[1],jointdata[2]))
                    linkid = jointdata[2]
                    deleteat!(jointlist,i)
                    break
                end
            end
            if linkid == origin.id
                break
            end
        end

        jointlist = [(joints[i].id,joints[i].parentid,joints[i].childids) for i=1:length(joints)]
        linkid = link2.id
        joint1id = 0
        joint2id = 0
        foundflag = false

        while true # check which predecessor link of link2 is also a predecessor link of link1
            for (i,jointdata) in enumerate(jointlist)
                if linkid ∈ jointdata[3]
                    joint2id = jointdata[1]
                    linkid = jointdata[2]
                    deleteat!(jointlist,i)
                    break
                end
            end
            for el in predlist
                if linkid == el[2]
                    joint1id = el[1]
                    foundflag = true
                    break
                end
            end
            foundflag && break            
        end

        # Find and remove joints to combine them
        joint1 = 0
        joint2 = 0
        for (i,joint) in enumerate(joints)
            if joint.id == joint1id
                joint1 = joint
                deleteat!(joints,i)
                break
            end
        end
        # if joint1id == joint2id # already combined joint
        #     joint2 = joint
        for (i,joint) in enumerate(joints)
            if joint.id == joint2id
                joint2 = joint
                deleteat!(joints,i)
                break
            end
        end
        
        joint = cat(joint1,joint2)
        push!(joints,joint)
        loopjoint = parse_loop_joint(xloopjoint, link1, link2, T)
        push!(loopjoints, loopjoint)
    end

    return joints, loopjoints
end 

function parse_urdf(filename, floating, ::Type{T}) where T
    xdoc = LightXML.parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xlinks = get_elements_by_tagname(xroot, "link")
    xjoints = get_elements_by_tagname(xroot, "joint")
    xloopjoints = get_elements_by_tagname(xroot, "loop_joint")

    materialdict = parse_robotmaterials(xroot, T)
    ldict = parse_links(xlinks, materialdict, T)
    origin, links, joints = parse_joints(xjoints, ldict, floating, T)

    joints, loopjoints = parse_loop_joints(xloopjoints, origin, joints, ldict, T)

    free(xdoc)

    return origin, links, joints, loopjoints
end


### After parsing

function set_parsed_values!(mechanism::Mechanism{T}, loopjoints) where T
    graph = mechanism.graph
    xjointlist = Dict{Int64,SVector{3,T}}() # stores id, x in world frame
    qjointlist = Dict{Int64,UnitQuaternion{T}}() # stores id, q in world frame

    for id in graph.rdfslist # from root to leaves
        component = getcomponent(mechanism, id)
        !(component isa Body) && continue # only for bodies

        body = component
        xbodylocal = body.state.xc
        qbodylocal = body.state.qc
        shape = body.shape

        parentid = get_parentid(mechanism, id, loopjoints)
        constraint = geteqconstraint(mechanism, parentid)

        grandparentid = constraint.parentid
        if grandparentid === nothing # predecessor is origin
            parentbody = mechanism.origin

            xparentbody = SA{T}[0; 0; 0]
            qparentbody = one(UnitQuaternion{T})

            xparentjoint = SA{T}[0; 0; 0]
            qparentjoint = one(UnitQuaternion{T})
        else
            parentbody = getbody(mechanism, grandparentid)

            grandgrandparentid = get_parentid(mechanism, grandparentid, loopjoints)
            parentconstraint = geteqconstraint(mechanism, grandgrandparentid)

            xparentbody = parentbody.state.xc # in world frame
            qparentbody = parentbody.state.qc # in world frame

            xparentjoint = xjointlist[parentconstraint.id] # in world frame
            qparentjoint = qjointlist[parentconstraint.id] # in world frame
        end

        ind1 = findfirst(x->x==id,constraint.childids)
        ind2 = ind1+1

        # urdf joint's x and q in parent's (parentbody) frame
        xjointlocal = vrotate(xparentjoint + vrotate(constraint.constraints[ind1].vertices[1], qparentjoint) - xparentbody, inv(qparentbody))
        qjointlocal = qparentbody \ qparentjoint * constraint.constraints[ind2].qoffset

        # store joint's x and q in world frame
        xjoint = xparentbody + vrotate(xjointlocal, qparentbody)
        qjoint = qparentbody * qjointlocal
        xjointlist[constraint.id] = xjoint
        qjointlist[constraint.id] = qjoint

        # difference to parent body (parentbody)
        qoffset = qjointlocal * qbodylocal

        # actual joint properties
        p1 = xjointlocal # in parent's (parentbody) frame
        p2 = vrotate(-xbodylocal, inv(qbodylocal)) # in body frame (xbodylocal and qbodylocal are both relative to the same (joint) frame -> rotationg by inv(body.q) gives body frame)
        constraint.constraints[ind1].vertices = (p1, p2)

        V3 = vrotate(constraint.constraints[ind2].V3', qjointlocal) # in parent's (parentbody) frame
        V12 = (svd(skew(V3)).Vt)[1:2,:]
        constraint.constraints[ind2].V3 = V3'
        constraint.constraints[ind2].V12 = V12
        constraint.constraints[ind2].qoffset = qoffset # in parent's (parentbody) frame

        # actual body properties
        setPosition!(body) # set everything to zero
        setPosition!(parentbody, body, p1 = p1, p2 = p2, Δq = qoffset)
        xbody = body.state.xc
        qbody = body.state.qc

        # shape relative
        if !(typeof(shape) <: EmptyShape)
            shape.xoffset = vrotate(xjoint + vrotate(shape.xoffset, qjoint) - xbody, inv(qbody))
            shape.qoffset = qoffset \ qjointlocal * shape.qoffset
        end
    end
    for (i,constraint) in enumerate(loopjoints)
        @assert length(constraint.childids) == 2 # Loop joint connects only two bodies

        parentid1 = constraint.parentid
        parentid2 = constraint.childids[1]
        if parentid1 === nothing # predecessor is origin
            parentbody1 = mechanism.origin

            xparentbody1 = SA{T}[0; 0; 0]
            qparentbody1 = one(UnitQuaternion{T})

            xparentjoint1 = SA{T}[0; 0; 0]
            qparentjoint1 = one(UnitQuaternion{T})
        else
            parentbody1 = getbody(mechanism, parentid1)

            grandparentid1 = get_parentid(mechanism, parentid1, loopjoints)
            parentconstraint1 = geteqconstraint(mechanism, grandparentid1)

            xparentbody1 = parentbody1.state.xc # in world frame
            qparentbody1 = parentbody1.state.qc # in world frame

            xparentjoint1 = xjointlist[parentconstraint1.id] # in world frame
            qparentjoint1 = qjointlist[parentconstraint1.id] # in world frame
        end
        parentbody2 = getbody(mechanism, parentid2)

        grandparentid2 = get_parentid(mechanism, parentid2, loopjoints)
        parentconstraint2 = geteqconstraint(mechanism, grandparentid2)

        xparentbody2 = parentbody2.state.xc # in world frame
        qparentbody2 = parentbody2.state.qc # in world frame

        xparentjoint2 = xjointlist[parentconstraint2.id] # in world frame
        qparentjoint2 = qjointlist[parentconstraint2.id] # in world frame


        ind1 = 1
        ind2 = ind1+1


        # urdf joint's x and q in parent's (parentbody) frame
        xjointlocal1 = vrotate(xparentjoint1 + vrotate(constraint.constraints[ind1].vertices[1], qparentjoint1) - xparentbody1, inv(qparentbody1))
        xjointlocal2 = vrotate(xparentjoint2 + vrotate(constraint.constraints[ind1].vertices[2], qparentjoint2) - xparentbody2, inv(qparentbody2))
        qjointlocal1 = qparentbody1 \ qparentjoint1 * constraint.constraints[ind2].qoffset

        # difference to parent body (parentbody)
        qoffset1 = qjointlocal1 * qparentbody2 #  qparentbody2 = body in for loop above

        # actual joint properties
        p1 = xjointlocal1 # in parent's (parentbody1) frame
        p2 = xjointlocal2 # in parent's (parentbody2) frame
        constraint.constraints[ind1].vertices = (p1, p2)

        V3 = vrotate(constraint.constraints[ind2].V3', qjointlocal1) # in parent's (parentbody1) frame
        V12 = (svd(skew(V3)).Vt)[1:2,:]
        constraint.constraints[ind2].V3 = V3'
        constraint.constraints[ind2].V12 = V12
        constraint.constraints[ind2].qoffset = qoffset1 # in parent's (parentbody1) frame
    end
end

function get_parentid(mechanism, id, loopjoints)
    graph = mechanism.graph
    conns = connections(graph, id)
    for connsid in conns
        constraint = geteqconstraint(mechanism, connsid)
        if constraint ∉ loopjoints && id ∈ constraint.childids
            return connsid
        end
    end

    return nothing
end