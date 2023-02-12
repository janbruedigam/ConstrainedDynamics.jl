import mujoco_py
from mujoco_py import load_model_from_path, MjSim, MjViewer
import os
from math import pi
import time
# Load the model from an XML file
model = load_model_from_path("n_snake_ball.xml")

# Create a simulation
sim = MjSim(model)
sim.model.opt.timestep = 0.01

# viewer = MjViewer(sim)

sim_state = sim.get_state()

N = 37
mindepth = 0.25
mintime = 100

for j in range(100):
    sim.set_state(sim_state)
    start_time = time.perf_counter()
    for i in range(500):
        sim.step()
        # for j in range(N-1):
        #     mindepth = min(mindepth,sim.data.body_xpos[j+1,2])
        # print(mindepth)
        # viewer.render()
    end_time = time.perf_counter()
    mintime = min(mintime,end_time - start_time)
    

print("Elapsed time: ", end_time - start_time)

# print(mindepth)