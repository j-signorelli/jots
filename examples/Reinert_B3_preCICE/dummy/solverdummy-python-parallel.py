#!/usr/bin/env python3
#from __future__ import division

# Updated from ajaust/precice-parallel-solverdummies

import numpy as np

from mpi4py import MPI

import precice

def Reinert_B3_Analytical(x,t,N=100):
    alpha = 2.5e-6
    L = 0.01
    k1 = 10
    k2 = 100
    T1 = 300
    T2 = 1300
    
    q_dot = 7.5e5
    T0 = 300
    
    D = ((k2-k1)/(T2-T1))*(1/(2*k1))
    theta_0 = (T0 - T1) + D*(T0-T1)**2
    
    theta = alpha*t/L**2 + 1/3 - x/L + 0.5*(x/L)**2
    for n in range(1,N+1):
        A = (2/np.pi**2)*(1/n**2)
        B = -n**2*np.pi**2*alpha*t/L**2
        C = n*np.pi*x/L
        theta -= A*np.exp(B)*np.cos(C)
        
    
    theta *= q_dot*L/k1
    theta += theta_0

    
    a = D
    b = 1-2*D*T1
    c = D*T1**2-T1-theta
    
    T = (-b+(b**2-4*a*c)**0.5)/(2*a)

    return T

def main():

    configuration_file_name = "../precice-config.xml"
    participant_name = "Dummy"

    write_data_name = 'Heat-Flux'
    read_data_name = 'Temperature'
    mesh_name = 'Dummy-Mesh'

    num_vertices = 1 # Only one needed at x=0 since problem is 1D

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
        vertices[i,0] = x[i]
        vertices[i,1] = 0
        vertices[i,2] = 0

    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)
    read_data_id = interface.get_data_id(read_data_name, mesh_id)
    write_data_id = interface.get_data_id(write_data_name, mesh_id)


    dt = interface.initialize()
    
    time = 0
    saved_time = 0
    write_data = np.array([-7.5e5]) # Write constant HF

    if (interface.is_action_required(precice.action_write_initial_data())):
        interface.write_block_scalar_data(write_data_id, vertex_ids, write_data)
        interface.mark_action_fulfilled(precice.action_write_initial_data())

    interface.initialize_data()

    while interface.is_coupling_ongoing():

        if interface.is_read_data_available():
            read_data = interface.read_block_scalar_data(read_data_id, vertex_ids)
        print("DUMMY ({}): Advancing in time".format(solver_process_index))
        
        time += dt


        T_received = read_data[0]
        T_analytical = Reinert_B3_Analytical(0, time)
        percent_error = abs(T_analytical - T_received)/T_analytical*100
        print("DUMMY ({}): t = {:.4f}, Temperature Percent Error is {}".format(solver_process_index, time, percent_error))

        if interface.is_write_data_required(dt):
            interface.write_block_scalar_data(
                write_data_id, vertex_ids, write_data)

        dt = interface.advance(dt)
        
    interface.finalize()
    print("DUMMY ({}): Closing python solver dummy...".format(solver_process_index))

# Can run only from terminal
if __name__ == '__main__':
    main()
