#pragma once
#include <sstream>
#include <cmath>

#include "mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "conductivity_model.hpp"

// Auxiliary class to prevent circular dependence, as PreciceAdapter is friend of PreciceBC
class PreciceAdapter;

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

class UniformConstantIsothermalBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformConstantIsothermalBC(const int attr, const double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::ISOTHERMAL, const_value){};
        bool IsEssential() const { return true; };
        std::string GetInitString() const;
};

class UniformConstantHeatFluxBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformConstantHeatFluxBC(const int attr, const double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::HEATFLUX, const_value){};
        bool IsEssential() const { return false; };
        std::string GetInitString() const;
};

class UniformSinusoidalIsothermalBC : public BoundaryCondition
{   
    private:
    protected:
        const double amplitude;
        const double ang_freq;
        const double phase;
        const double vert_shift;

        const double& time_ref; // reference to the time

    public:
        UniformSinusoidalIsothermalBC(const int attr, const double& in_tref, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert) : BoundaryCondition(attr, BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL), time_ref(in_tref), amplitude(in_amp), ang_freq(in_angfreq), phase(in_phase), vert_shift(in_vert) {};
        bool IsConstant() const { return false; }
        bool IsEssential() const { return true; }
        void InitCoefficient();
        void UpdateCoeff();
        std::string GetInitString() const;

};

class PreciceBC : public BoundaryCondition
{
    friend class PreciceAdapter; // Allow PreciceAdapter access to everything

    private:
    protected:

        //const std::string TEMPERATURE = "Temperature";
        //const std::string HEATFLUX = "Heat-Flux";

        //const ConductivityModel* cond_model; // Ptr to conductivity model - not allocated here
        //const SolverState* curr_state; // Ptr to solver state -- not allocated here

        mfem::ParFiniteElementSpace& fespace;

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

        double* read_data_arr;
        double* write_data_arr;

        mfem::Array<double> coeff_dof_values;
        mfem::ParGridFunction* coeff_gf;

        bool update_flag;

        // Set by adapter:
        int mesh_id;
        int* vertex_ids;
        int read_data_id;
        int write_data_id;
    

        static void GetBdrTemperatures(const mfem::ParGridFunction* T_gf, const mfem::Array<int> in_bdr_elem_indices, double* nodal_temperatures); // Precondition: in_bdr_elem_indices length is correct
        static void GetBdrWallHeatFlux(const mfem::ParGridFunction* T_gf, const ConductivityModel* in_cond, const mfem::Array<int> in_bdr_elem_indices, double* nodal_wall_heatfluxes); // Precondition same above
        mutable mfem::ParGridFunction* temp_gf;
    public:
        PreciceBC(const int attr, const BOUNDARY_CONDITION in_type, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const bool is_restart, const double in_value, const std::string in_read, const std::string in_write);
        bool IsConstant() const { return false; };
        void InitCoefficient();
        void UpdateCoeff();

        virtual void RetrieveInitialWriteData(const mfem::Vector T, const ConductivityModel* cond_model) = 0;
        virtual void RetrieveWriteData(const mfem::Vector T, const ConductivityModel* cond_model) = 0;


        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
        ~PreciceBC();
       
};


class PreciceIsothermalBC : public PreciceBC
{
    private:

    protected:
    public:
        PreciceIsothermalBC(const int attr, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const bool is_restart, const double in_value) : PreciceBC(attr, BOUNDARY_CONDITION::PRECICE_ISOTHERMAL, f, in_mesh, is_restart, in_value, "Temperature", "Heat-Flux") {};
        void RetrieveInitialWriteData(const mfem::Vector T, const ConductivityModel* cond_model);
        void RetrieveWriteData(const mfem::Vector T, const ConductivityModel* cond_model);
        std::string GetInitString() const;

        bool IsEssential() const { return true; };
};

class PreciceHeatFluxBC : public PreciceBC
{
    private:
    
    protected:

    public:
        PreciceHeatFluxBC(const int attr, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const bool is_restart, const double in_value) : PreciceBC(attr, BOUNDARY_CONDITION::PRECICE_HEATFLUX, f, in_mesh, is_restart, in_value, "Heat-Flux", "Temperature") {};
        void RetrieveInitialWriteData(const mfem::Vector T, const ConductivityModel* cond_model);
        void RetrieveWriteData(const mfem::Vector T, const ConductivityModel* cond_model);
        std::string GetInitString() const;

        bool IsEssential() const { return false; };
};