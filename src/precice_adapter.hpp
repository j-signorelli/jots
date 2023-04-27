#pragma once

#include "precice/SolverInterface.hpp"

#include "solver_state.hpp"
#include "boundary_condition.hpp"

class PreciceAdapter
{
    private:
        precice::SolverInterface* interface;
        SolverState* saved_state;
        
        PreciceBC** precice_bcs; // Individual BC ptrs not allocated here
        size_t num_bcs;
        
        const std::string participant_name;
        const std::string config_file;
        const int rank;
        const int size;
        int dim;

    public:
        PreciceAdapter(const std::string in_part_name, const std::string in_config, const int r, const int s);

        void AddPreciceBCs(BoundaryCondition** in_bcs, std::vector<int> precice_bc_indices);
        
        double Initialize() { return interface->initialize(); };
        
        void SendInitialData();

        void InitializeData() { interface->initializeData(); };

        void SaveOldState(const SolverState* state);
        
        void ReloadOldState(SolverState* state) const;

        int GetDimension() const { return dim; };

        ~PreciceAdapter();



}