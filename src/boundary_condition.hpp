#pragma once
#include <sstream>

#include "mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "conductivity_model.hpp"

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

        precice::SolverInterface* interface; // Not allocated here

        mfem::ParFiniteElementSpace* fespace; // Not allocated here
        const mfem::ParGridFunction* T_gf; // Ptr to main solution GridFunction - Not allocated here, required for obtaining write data
        const ConductivityModel* cond_model; // Ptr to conductivity model - not allocated here

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
        
        double* readDataArr;
        double* writeDataArr;


        mfem::ParGridFunction coeff_gf;
        mfem::Array<double> coeff_values;

        static void GetBdrTemperatures(const mfem::ParGridFunction* T_gf, const mfem::Array<int> in_bdr_elem_indices, double* nodal_temperatures); // Precondition: in_bdr_elem_indices length is correct
        static void GetBdrWallHeatFlux(const mfem::ParGridFunction* T_gf, const ConductivityModel* cond_model, const mfem::Array<int> in_bdr_elem_indices, double* nodal_wall_heatfluxes); // Precondition same above

    public:
        preCICEBC(int attr, BOUNDARY_CONDITION in_type, precice::SolverInterface* in, const mfem::ParGridFunction* in_T_gf, ConductivityModel* in_cond, bool is_restart, std::string mesh_name, double initial_value, std::string read_data_name, std::string write_data_name);
        bool IsConstant() const { return false; };
        
        void InitCoefficient();
        void UpdateCoeff();

        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
        ~preCICEBC();
       
};


class preCICEIsothermalBC : public preCICEBC
{
    private:

    protected:

    public:
        preCICEIsothermalBC(int attr, precice::SolverInterface* in, mfem::ParGridFunction* in_T_gf, ConductivityModel* in_cond, bool is_restart, std::string mesh_name, double initial_value) : preCICEBC(attr, BOUNDARY_CONDITION::PRECICE_ISOTHERMAL, in, in_T_gf, in_cond, is_restart, mesh_name, initial_value, "Temperature", "Heat-Flux") {}; // TODO: Generalize strings!!!
        bool IsEssential() const { return true; };
        std::string GetInitString() const {return "";};
};

class preCICEHeatFluxBC : public preCICEBC
{
    private:
    
    protected:

    public:
        preCICEHeatFluxBC(int attr, precice::SolverInterface* in, mfem::ParGridFunction* in_T_gf, ConductivityModel* in_cond, bool is_restart, std::string mesh_name, double initial_value) : preCICEBC(attr, BOUNDARY_CONDITION::PRECICE_HEATFLUX, in, in_T_gf, in_cond, is_restart, mesh_name, initial_value, "Heat-Flux", "Temperature") {};
        bool IsEssential() const { return false; };
        std::string GetInitString() const {return "";};
};