import mujoco_py
from mujoco_py import load_model_from_path, MjSim, MjViewer
import os
from math import pi
import time
import numpy

# Load the model from an XML file
model = load_model_from_path("box_toss.xml")

# Create a simulation
sim = MjSim(model)
sim.model.opt.timestep = 0.01

viewer = MjViewer(sim)

# Run the simulation for some number of steps
joint_id1 = model.joint_name2id("joint1")

sim_state = sim.get_state()


for d in [0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25]:
#  for d in [0.65]:
    sim.set_state(sim_state)
    sim.data.qpos[joint_id1:joint_id1+3] = [0,0,d]
    depth = numpy.zeros(100)

    for i in range(100):
        depth[i] = sim.data.qpos[joint_id1+2]
        # print(depth[i])
        sim.step()
        # viewer.render()

    # if os.getenv('TESTING') is not None:
    #     break
    
    print(min(depth))