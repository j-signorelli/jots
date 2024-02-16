#pragma once

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "boundary_condition.hpp"
#include "material_property.hpp"

class PreciceAdapter
{
    private:
        precice::SolverInterface interface;
        
        PreciceBC** precice_bcs; // Individual BC ptrs not allocated here
        size_t num_bcs;
        
        int dim;

        mfem::Vector old_state_u;

    public:

        static const std::string cowid;
        static const std::string cowic;
        static const std::string coric;

        PreciceAdapter(const std::string in_part_name, const std::string in_config, const int r, const int s, MPI_Comm comm=MPI_COMM_WORLD);

        void SetPreciceBCs(BoundaryCondition** in_bcs, std::vector<int> precice_bc_indices);
        
        precice::SolverInterface& Interface() { return interface; };

        void GetReadData();

        void WriteData(const mfem::ParGridFunction &u_gf);

        void SaveOldState(const mfem::Vector u);

        void ReloadOldState(mfem::Vector& u) const;

        int GetDimension() const { return dim; };

        ~PreciceAdapter();



};