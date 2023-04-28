#!/usr/bin/env python3
#from __future__ import division

# Updated from ajaust/precice-parallel-solverdummies

import numpy as np

from mpi4py import MPI

import precice

configuration_file_name = "../precice-config.xml"
participant_name = "Dummy"

write_data_name = 'Temperature'
read_data_name = 'Heat-Flux'
mesh_name = 'Dummy-Mesh'

num_vertices = 100  # Number of vertices

comm = MPI.COMM_WORLD
solver_process_index = comm.Get_rank()
solver_process_size = comm.Get_size()

interface = precice.Interface(participant_name, configuration_file_name,
                              solver_process_index, solver_process_size)

mesh_id = interface.get_mesh_id(mesh_name)
dimensions = interface.get_dimensions()

vertices = np.zeros((num_vertices, dimensions))
read_data = np.zeros((num_vertices))
write_data = np.zeros((num_vertices))

x = np.linspace(0,1,num_vertices)
for i in range(num_vertices):
    vertices[i,0] = x[i]
    vertices[i,1] = -0.25
    read_data[i] = i+solver_process_index
    write_data[i] = 310

vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)
read_data_id = interface.get_data_id(read_data_name, mesh_id)
write_data_id = interface.get_data_id(write_data_name, mesh_id)


dt = interface.initialize()

if (interface.is_action_required(precice.action_write_initial_data())):
    interface.write_block_scalar_data(write_data_id, vertex_ids, write_data)
    interface.mark_action_fulfilled(precice.action_write_initial_data())

interface.initialize_data()

while interface.is_coupling_ongoing():
    if interface.is_action_required(
            precice.action_write_iteration_checkpoint()):
        print("DUMMY ({}): Writing iteration checkpoint".format(solver_process_index))
        interface.mark_action_fulfilled(
            precice.action_write_iteration_checkpoint())

    if interface.is_read_data_available():
        read_data = interface.read_block_scalar_data(read_data_id, vertex_ids)


    if interface.is_write_data_required(dt):
        interface.write_block_scalar_data(
            write_data_id, vertex_ids, write_data)

    print("DUMMY ({}): Advancing in time".format(solver_process_index))
    dt = interface.advance(dt)

    if interface.is_action_required(
            precice.action_read_iteration_checkpoint()):
        print("DUMMY ({}): Reading iteration checkpoint".format(solver_process_index))
        interface.mark_action_fulfilled(
            precice.action_read_iteration_checkpoint())

interface.finalize()
print("DUMMY ({}): Closing python solver dummy...".format(solver_process_index))