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
        int bdr_attr;
        BOUNDARY_CONDITION bc_type;
        mfem::Coefficient* coeff;
        
    public:
        BoundaryCondition(int attr, BOUNDARY_CONDITION in_type) : bdr_attr(attr), bc_type(in_type) {};
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
        double uniform_value;
    public:
        UniformConstantBC(int attr, BOUNDARY_CONDITION in_type, double in_value) : BoundaryCondition(attr, in_type), uniform_value(in_value) {};
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
        UniformIsothermalBC(int attr, double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::ISOTHERMAL, const_value){};
        bool IsEssential() const { return true; };
        std::string GetInitString() const;
};

class UniformHeatFluxBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformHeatFluxBC(int attr, double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::HEATFLUX, const_value){};
        bool IsEssential() const { return false; };
        std::string GetInitString() const;
};


class preCICEBC : public BoundaryCondition
{
    private:
    protected:

        const std::string TEMPERATURE = "Temperature";
        const std::string HEATFLUX = "Heat-Flux";

        precice::SolverInterface* interface; // Not allocated here

        const ConductivityModel* cond_model; // Ptr to conductivity model - not allocated here
        const SolverState* curr_state; // Ptr to solver state -- not allocated here

        std::string readDataName;
        std::string writeDataName;

        int meshID;
        double* coords;
        int* vertexIDs;
        int readDataID;
        int writeDataID;

        bool restart;

        mfem::Array<int> bdr_elem_indices;
        mfem::Array<int> bdr_dof_indices; // Includes non-true DOFs
        int num_dofs;

        double* readDataArr;
        double* writeDataArr;

        double initial_value;

        mfem::ParGridFunction* coeff_gf;
        mfem::Array<double> coeff_values;

        static void GetBdrTemperatures(const mfem::ParGridFunction* T_gf, const mfem::Array<int> in_bdr_elem_indices, double* nodal_temperatures); // Precondition: in_bdr_elem_indices length is correct
        static void GetBdrWallHeatFlux(const mfem::ParGridFunction* T_gf, const ConductivityModel* in_cond, const mfem::Array<int> in_bdr_elem_indices, double* nodal_wall_heatfluxes); // Precondition same above

    public:
        preCICEBC(int attr, BOUNDARY_CONDITION in_type, const SolverState* in_state, precice::SolverInterface* in_interface, const ConductivityModel* in_cond, bool is_restart, std::string mesh_name, double in_value, std::string read_data_name, std::string write_data_name);
        bool IsConstant() const { return false; };
        void InitCoefficient();
        void UpdateCoeff();

        virtual void GetInitialWriteDataFxn() = 0; // Function to get initial read data from restarted state of T_gf
        virtual void GetWriteDataFxn() = 0;


        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
        ~preCICEBC();
       
};


class preCICEIsothermalBC : public preCICEBC
{
    private:

    protected:

    public:
        preCICEIsothermalBC(int attr, const SolverState* in_state, precice::SolverInterface* in_interface, const ConductivityModel* in_cond, bool is_restart, std::string mesh_name, double in_value) : preCICEBC(attr, BOUNDARY_CONDITION::PRECICE_ISOTHERMAL, in_state, in_interface, in_cond, is_restart, mesh_name, in_value, "Temperature", "Heat-Flux") {}; // TODO: Understand why I cannot do as done below. I get bad_alloc error
        void GetInitialWriteDataFxn() { GetBdrTemperatures(curr_state->GetGF(), bdr_elem_indices, readDataArr); }; // SHOULD BE HEAT FLUX
        void GetWriteDataFxn() { GetBdrWallHeatFlux(curr_state->GetGF(), cond_model, bdr_elem_indices, writeDataArr); };
        bool IsEssential() const { return true; };
        std::string GetInitString() const {return "";};
};

class preCICEHeatFluxBC : public preCICEBC
{
    private:
    
    protected:

    public:
        preCICEHeatFluxBC(int attr, const SolverState* in_state, precice::SolverInterface* in_interface, const ConductivityModel* in_cond, bool is_restart, std::string mesh_name, double in_value) : preCICEBC(attr, BOUNDARY_CONDITION::PRECICE_HEATFLUX, in_state, in_interface, in_cond, is_restart, mesh_name, in_value, HEATFLUX, TEMPERATURE) {};
        void GetInitialWriteDataFxn() { GetBdrWallHeatFlux(curr_state->GetGF(), cond_model, bdr_elem_indices, readDataArr); }; // SHOULD BE TEMPERATURE
        void GetWriteDataFxn() { GetBdrTemperatures(curr_state->GetGF(), bdr_elem_indices, readDataArr); };
        
        bool IsEssential() const { return false; };
        std::string GetInitString() const {return "";};
};