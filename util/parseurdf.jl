# Some adaptions from RigidBodyDynamics

function parse_scalar(::Type{T}, e::XMLElement, name::String) where {T}
    parse(T, attribute(e, name))
end

function parse_scalar(::Type{T}, e::XMLElement, name::String, default::String) where {T}
    parse(T, e == nothing ? default : attribute(e, name))
end

function parse_vector(::Type{T}, e::Union{XMLElement, Nothing}, name::String, default::String) where {T}
    usedefault = e == nothing || attribute(e, name) == nothing # TODO: better handling of required attributes
    [parse(T, str) for str in split(usedefault ? default : attribute(e, name))]
end

function parse_inertia_matrix(::Type{T}, xml_inertia::Nothing) where T
    zeros(T,3,3)
end

function parse_inertia_matrix(::Type{T}, xml_inertia::XMLElement) where T
    ixx = parse_scalar(T, xml_inertia, "ixx", "0")
    ixy = parse_scalar(T, xml_inertia, "ixy", "0")
    ixz = parse_scalar(T, xml_inertia, "ixz", "0")
    iyy = parse_scalar(T, xml_inertia, "iyy", "0")
    iyz = parse_scalar(T, xml_inertia, "iyz", "0")
    izz = parse_scalar(T, xml_inertia, "izz", "0")
    [ixx ixy ixz; ixy iyy iyz; ixz iyz izz]
end

function parse_pose(::Type{T}, xml_pose::Nothing) where T
    x = zero(T,3)
    q = Quaternion{T}()
    x, q
end

function parse_pose(::Type{T}, xml_pose::XMLElement) where T
    x = parse_vector(T, xml_pose, "xyz", "0 0 0")
    rpy = parse_vector(T, xml_pose, "rpy", "0 0 0")
    q = Quaternion(RotZYX(rpy[3], rpy[2], rpy[1]))
    x, q
end

function find_parent_joint(link::XMLLink, joints::Vector{XMLJoint})
    parentid = 0
    for (i,joint) in enumerate(joints)
        if link.name == joint.child
            parentid = i
        end
    end
    return parentid
end

function find_child_joint(link::XMLLink, joints::Vector{XMLJoint})
    childid = Vector{Int64}(undef,0)
    for (i,joint) in enumerate(joints)
        if link.name == joint.parent
            push!(childid,i)
        end
    end
    return childid
end

function parse_joint_type(::Type{T}, xml_joint::XMLElement) where T
    urdf_joint_type = attribute(xml_joint, "type")
    if urdf_joint_type == "revolute" || urdf_joint_type == "continuous" || urdf_joint_type == "prismatic"
        axis = parse_vector(T, find_element(xml_joint, "axis"), "xyz", "1 0 0")
        dof = 1
        return "revolute", dof, axis
    # elseif urdf_joint_type == "floating" || urdf_joint_type == "fixed"
    #     return joint_type{T}()
    # elseif urdf_joint_type == "planar"
    #     urdf_axis = SVector{3}(parse_vector(T, find_element(xml_joint, "axis"), "xyz", "1 0 0"))
    #     # The URDF spec says that a planar joint allows motion in a
    #     # plane perpendicular to the axis.
    #     R = Rotations.rotation_between(SVector(0, 0, 1), urdf_axis)
    #     x_axis = R * SVector(1, 0, 0)
    #     y_axis = R * SVector(0, 1, 0)
    #     return joint_type(x_axis, y_axis)
    else
        error("joint type $(urdf_joint_type) not recognized")
    end
end

function parse_inertia(::Type{T}, xml_inertial::XMLElement) where T
    x, q = parse_pose(T, find_element(xml_inertial, "origin"))
    inertia = parse_inertia_matrix(T, find_element(xml_inertial, "inertia"))
    mass = parse_scalar(T, find_element(xml_inertial, "mass"), "value", "0")

    x, q, mass, inertia
end

function parse_connected_links(xml_joint::XMLElement)
    xml_parent = find_element(xml_joint, "parent")
    parent = attribute(xml_parent, "link")
    xml_child = find_element(xml_joint, "child")
    child = attribute(xml_child, "link")
    parent, child
end

function parse_link(::Type{T}, xml_link::XMLElement) where T
    xml_inertial = find_element(xml_link, "inertial")
    if xml_inertial!=nothing
        x, q, mass, inertia = parse_inertia(T, xml_inertial)
    else
        x = zeros(T,3)
        q = Quaternion{T}()
        mass = zero(T)
        inertia = zeros(T,3,3)
    end
    name = attribute(xml_link, "name")
    XMLLink{T}(name,0,x,q,mass,inertia,0,[-x])
end

function parse_joint(::Type{T}, xml_joint::XMLElement) where T
    name = attribute(xml_joint, "name")
    type, dof, axis = parse_joint_type(T, xml_joint)
    x, q = parse_pose(T, find_element(xml_joint, "origin"))
    parent, child = parse_connected_links(xml_joint)

    XMLJoint{T}(name,type,x,q,axis,parent,child,0,0,(0,1),dof)
end


function parse_urdf(filename::AbstractString; scalar_type::Type{T}=Float64,
    tend::T=10., dt::T=.01, g::T=-9.81, rootid=0) where T
    xdoc = parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xml_links = get_elements_by_tagname(xroot, "link")
    xml_joints = get_elements_by_tagname(xroot, "joint")

    linkscoll = Vector{XMLLink}(undef,0)
    jointscoll = Vector{XMLJoint}(undef,0)

    for xml_link in xml_links
        link = parse_link(T, xml_link)
        push!(linkscoll,link)
    end
    for xml_joint in xml_joints
        push!(jointscoll,parse_joint(T, xml_joint))
    end

    # free(xdoc)

    links = Vector{Link{T}}(undef,0)
    constraints = Vector{Constraint{T}}(undef,0)
    origin = 0
    originid = 0

    for (i,link) in enumerate(linkscoll)
        childid = find_child_joint(link, jointscoll)
        for ind in childid
            joint = jointscoll[ind]
            push!(link.p,link.p[1]+joint.x)
            joint.pids = (length(link.p),joint.pids[2])
            joint.parentid = i
        end

        parentid = find_parent_joint(link, jointscoll)
        if parentid == 0
            link.dof = 0
            originid = i
            origin = Link{T}(link)
        else
            joint = jointscoll[parentid]
            link.dof = joint.dof
            joint.childid = i

            push!(links,Link{T}(link))
        end
    end


    for joint in jointscoll
        parentid = joint.parentid
        if parentid<originid
            link1=links[parentid]
        elseif parentid==originid
            link1=origin
        elseif parentid>=originid
            link1=links[parentid-1]
        end
        childid = joint.childid
        if childid<originid
            link2 = links[childid]
        elseif childid==originid
            error("sth is wrong")
        elseif childid>=originid
            link2 = links[childid-1]
        end

        push!(constraints,Combined(joint,link1,link2))
    end

    Robot(origin,links,constraints,tend=tend, dt=dt, g=g, rootid=rootid)
end
