#!/usr/bin/env python3
#from __future__ import division

# Updated from ajaust/precice-parallel-solverdummies

import numpy as np

from mpi4py import MPI

import precice

def Reinert_B1_Analytical_HF(x,t,N=100):
    alpha = 2.5e-6
    L = 0.01
    k = 10
    dTdx = 0
    T_0 = 300
    T_D = 500

    for n in range(0,N+1):
        A = (-1)**n/(2*n+1)
        B = (-(2*n+1)**2*np.pi**2*alpha*t)/(4*L**2)
        C = (2*n+1)*np.pi/(2*L)
        dTdx = dTdx + (T_D - T_0)*C*(4/np.pi)*A*np.exp(B)*np.sin(C*x)
    
    return -k*dTdx

def main():

    configuration_file_name = "../precice-config.xml"
    participant_name = "Dummy"

    write_data_name = 'Temperature'
    read_data_name = 'Heat-Flux'
    mesh_name = 'Dummy-Mesh'

    num_vertices = 1 # Only one needed at x=0.01 since problem is 1D

    comm = MPI.COMM_WORLD
    solver_process_index = comm.Get_rank()
    solver_process_size = comm.Get_size()

    interface = precice.Interface(participant_name, configuration_file_name,
                                solver_process_index, solver_process_size)

    mesh_id = interface.get_mesh_id(mesh_name)
    dimensions = interface.get_dimensions() #3D in actuality

    vertices = np.zeros((num_vertices, dimensions))
    read_data = np.zeros((num_vertices)) # this is always ignored
    write_data = np.zeros((num_vertices))
    
    x = np.linspace(0,1,num_vertices)
    for i in range(num_vertices):
        vertices[i,0] = 0.01 #x[i]
        vertices[i,1] = 0
        vertices[i,2] = 0

    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)
    read_data_id = interface.get_data_id(read_data_name, mesh_id)
    write_data_id = interface.get_data_id(write_data_name, mesh_id)


    dt = interface.initialize()
    
    time = 0
    write_data = np.array([500]) # Write constant temperature

    if (interface.is_action_required(precice.action_write_initial_data())):
        interface.write_block_scalar_data(write_data_id, vertex_ids, write_data)
        interface.mark_action_fulfilled(precice.action_write_initial_data())

    interface.initialize_data()

    while interface.is_coupling_ongoing():

        if interface.is_read_data_available():
            read_data = interface.read_block_scalar_data(read_data_id, vertex_ids)
        print("DUMMY ({}): Advancing in time".format(solver_process_index))
        
        time += dt


        q_received = read_data[0]
        q_analytical = Reinert_B1_Analytical_HF(0.01, time)
        percent_error = abs((q_analytical - q_received)/q_analytical)*100
        print("DUMMY ({}): t = {:.4f}, HF Percent Error is {}".format(solver_process_index, time, percent_error))

        if interface.is_write_data_required(dt):
            interface.write_block_scalar_data(
                write_data_id, vertex_ids, write_data)

        dt = interface.advance(dt)
        
    interface.finalize()
    print("DUMMY ({}): Closing python solver dummy...".format(solver_process_index))

# Can run only from terminal
if __name__ == '__main__':
    main()
