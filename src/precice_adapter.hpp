#pragma once

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "boundary_condition.hpp"
#include "conductivity_model.hpp"

class PreciceAdapter
{
    private:
        precice::SolverInterface* interface;
        
        PreciceBC** precice_bcs; // Individual BC ptrs not allocated here
        size_t num_bcs;
        
        const std::string participant_name;
        const std::string config_file;
        const int rank;
        const int size;
        int dim;

        mfem::Vector old_state_T;

    public:

        static const std::string cowid;
        static const std::string cowic;
        static const std::string coric;

        PreciceAdapter(const std::string in_part_name, const std::string in_config, const int r, const int s);

        void AddPreciceBCs(BoundaryCondition** in_bcs, std::vector<int> precice_bc_indices);
        
        precice::SolverInterface* Interface() { return interface; };

        void WriteInitialData(const mfem::Vector T, const ConductivityModel* cond_model);

        void GetReadData();

        void WriteData(const mfem::Vector T, const ConductivityModel* cond_model);

        void SaveOldState(const mfem::Vector T);

        void ReloadOldState(mfem::Vector& T) const;

        int GetDimension() const { return dim; };

        ~PreciceAdapter();



};