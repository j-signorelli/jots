#pragma once
#include <sstream>

#include "mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "conductivity_model.hpp"
#include "solver_state.hpp"

class BoundaryCondition
{
    private:
    protected:
        const int bdr_attr;
        const BOUNDARY_CONDITION bc_type;
        mfem::Coefficient* coeff;
        
    public:
        BoundaryCondition(const int attr, const BOUNDARY_CONDITION in_type) : bdr_attr(attr), bc_type(in_type) {};
        int GetBdrAttr() const { return bdr_attr; }
        BOUNDARY_CONDITION GetType() const { return bc_type; };
        mfem::Coefficient* GetCoeffPtr() const { return coeff; };

        virtual void InitCoefficient() = 0;
        virtual void UpdateCoeff() = 0;
        virtual bool IsEssential() const = 0;
        virtual bool IsConstant() const = 0; // true if d/dt is 0 for this coefficient
        virtual std::string GetInitString() const = 0;

        ~BoundaryCondition() {delete coeff;};
};

class UniformConstantBC : public BoundaryCondition
{   
    private:
    protected:
        const double uniform_value;
    public:
        UniformConstantBC(const int attr, const BOUNDARY_CONDITION in_type, const double in_value) : BoundaryCondition(attr, in_type), uniform_value(in_value) {};
        bool IsConstant() const { return true; };
        double GetValue() const { return uniform_value; };
        void InitCoefficient() { coeff = new mfem::ConstantCoefficient(uniform_value); };
        void UpdateCoeff() {};
        virtual std::string GetInitString() const = 0;
};

class UniformIsothermalBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformIsothermalBC(const int attr, const double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::ISOTHERMAL, const_value){};
        bool IsEssential() const { return true; };
        std::string GetInitString() const;
};

class UniformHeatFluxBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformHeatFluxBC(const int attr, const double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::HEATFLUX, const_value){};
        bool IsEssential() const { return false; };
        std::string GetInitString() const;
};


class PreciceBC : public BoundaryCondition
{
    private:
    protected:

        //const std::string TEMPERATURE = "Temperature";
        //const std::string HEATFLUX = "Heat-Flux";

        //const ConductivityModel* cond_model; // Ptr to conductivity model - not allocated here
        //const SolverState* curr_state; // Ptr to solver state -- not allocated here

        const mfem::ParFiniteElementSpace& fespace;
        const std::string mesh_name;
        const bool restart;
        const double default_value;
        const std::string read_data_name;
        const std::string write_data_name;
        const int dim;
        
        mfem::Array<int> bdr_elem_indices;
        mfem::Array<int> bdr_dof_indices; // Includes non-true DOFs
        double* coords;
        int num_dofs;
        int* vertex_ids;

        double* read_data_arr;
        double* write_data_arr;

        mfem::Array<double> coeff_dof_values;
        mfem::ParGridFunction* coeff_gf;

        // Set by adapter:
        int mesh_id;
        int* vertex_ids;
        int read_data_id;
        int write_data_id;

        static void GetBdrTemperatures(const mfem::ParGridFunction* T_gf, const mfem::Array<int> in_bdr_elem_indices, double* nodal_temperatures); // Precondition: in_bdr_elem_indices length is correct
        static void GetBdrWallHeatFlux(const mfem::ParGridFunction* T_gf, const ConductivityModel* in_cond, const mfem::Array<int> in_bdr_elem_indices, double* nodal_wall_heatfluxes); // Precondition same above

    public:
        PreciceBC(const int attr, const BOUNDARY_CONDITION in_type, const mfem::ParFiniteElementSpace& f, const std::string in_mesh, const bool is_restart, const double in_value, const std::string in_read, const std::string in_write);
        bool IsConstant() const { return false; };
        void InitCoefficient();
        void UpdateCoeff();

        std::string GetMeshName() const { return mesh_name; };

        void SetMeshID(int id) { mesh_id = id; };
        int GetMeshID() const { return mesh_id; };

        double* GetCoords() const { return coords; };

        std::string GetReadDataName() const { return read_data_name; };
        std::string GetWriteDataName() const { return write_data_name; };
        
        void SetReadDataID(int id) { read_data_id = id; };
        void SetWriteDataID(int id) { write_data_id = id; };

        virtual void GetInitialWriteDataFxn() = 0; // Function to get initial read data from restarted state of T_gf
        virtual void GetWriteDataFxn() = 0;


        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
        ~preCICEBC();
       
};


class PreciceIsothermalBC : public PreciceBC
{
    private:

    protected:

    public:
        PreciceIsothermalBC(const int attr, const mfem::ParFiniteElementSpace& f, const std::string in_mesh, const bool is_restart, const double in_value) : PreciceBC(attr, BOUNDARY_CONDITION::PRECICE_ISOTHERMAL, f, in_mesh, is_restart, in_value, "Temperature", "Heat-Flux") {};
        void GetInitialWriteDataFxn() { GetBdrTemperatures(curr_state->GetGF(), bdr_elem_indices, readDataArr); }; // SHOULD BE HEAT FLUX
        void GetWriteDataFxn() { GetBdrWallHeatFlux(curr_state->GetGF(), cond_model, bdr_elem_indices, writeDataArr); };
        bool IsEssential() const { return true; };
};

class PreciceHeatFluxBC : public PreciceBC
{
    private:
    
    protected:

    public:
        PreciceHeatFluxBC(const int attr, const mfem::ParFiniteElementSpace& f, const std::string in_mesh, const bool is_restart, const double in_value) : PreciceBC(attr, BOUNDARY_CONDITION::PRECICE_HEATFLUX, f, in_mesh, is_restart, in_value, "Heat-Flux", "Temperature") {};
        void GetInitialWriteDataFxn() { GetBdrWallHeatFlux(curr_state->GetGF(), cond_model, bdr_elem_indices, readDataArr); }; // SHOULD BE TEMPERATURE
        void GetWriteDataFxn() { GetBdrTemperatures(curr_state->GetGF(), bdr_elem_indices, readDataArr); };
        
        bool IsEssential() const { return false; };
};