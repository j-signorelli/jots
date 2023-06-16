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
        mfem::Coefficient& GetCoeffRef() const { return *coeff; }; // To be used only for assigning to linear/bilinear forms or projecting coeff, so declared const

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
        UniformConstantBC(const int attr, const BOUNDARY_CONDITION in_type, const double in_value);
        bool IsConstant() const { return true; };
        double GetValue() const { return uniform_value; };
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

class UniformSinusoidalBC : public BoundaryCondition
{   
    private:
    protected:
        const double amplitude;
        const double ang_freq;
        const double phase;
        const double vert_shift;

        const double& time_ref; // reference to the time

    public:
        UniformSinusoidalBC(const int attr, BOUNDARY_CONDITION in_type, const double& in_tref, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert);
        bool IsConstant() const { return false; }
        void UpdateCoeff();
        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;

};

class UniformSinusoidalIsothermalBC : public UniformSinusoidalBC
{   
    private:
    protected:
    public:
        UniformSinusoidalIsothermalBC(const int attr, const double& in_tref, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert) : UniformSinusoidalBC(attr, BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL, in_tref, in_amp, in_angfreq, in_phase, in_vert) {};
        bool IsEssential() const { return true; }
        std::string GetInitString() const;

};

class UniformSinusoidalHeatFluxBC : public UniformSinusoidalBC
{   
    private:
    protected:
    public:
        UniformSinusoidalHeatFluxBC(const int attr, const double& in_tref, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert) : UniformSinusoidalBC(attr, BOUNDARY_CONDITION::SINUSOIDAL_HEATFLUX, in_tref, in_amp, in_angfreq, in_phase, in_vert) {};
        bool IsEssential() const { return false; }
        std::string GetInitString() const;

};

class PreciceBC : public BoundaryCondition
{
    friend class PreciceAdapter; // Allow PreciceAdapter access to everything

    private:
    protected:

        //const std::string TEMPERATURE = "Temperature";
        //const std::string HEATFLUX = "Heat-Flux";

        mfem::ParFiniteElementSpace& fespace;

        const std::string mesh_name;
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
        PreciceBC(const int attr, const BOUNDARY_CONDITION in_type, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const double in_value, const std::string in_read, const std::string in_write);
        bool IsConstant() const { return false; };
        void UpdateCoeff();

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
        PreciceIsothermalBC(const int attr, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const double in_value) : PreciceBC(attr, BOUNDARY_CONDITION::PRECICE_ISOTHERMAL, f, in_mesh, in_value, "Temperature", "Heat-Flux") {};
        void RetrieveWriteData(const mfem::Vector T, const ConductivityModel* cond_model);
        std::string GetInitString() const;

        bool IsEssential() const { return true; };
};

class PreciceHeatFluxBC : public PreciceBC
{
    private:
    
    protected:

    public:
        PreciceHeatFluxBC(const int attr, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const double in_value) : PreciceBC(attr, BOUNDARY_CONDITION::PRECICE_HEATFLUX, f, in_mesh, in_value, "Heat-Flux", "Temperature") {};
        void RetrieveWriteData(const mfem::Vector T, const ConductivityModel* cond_model);
        std::string GetInitString() const;

        bool IsEssential() const { return false; };
};