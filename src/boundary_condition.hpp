#pragma once
#include <sstream>
#include <cmath>

#include "mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "material_property.hpp"

// Auxiliary class to prevent circular dependence, as PreciceAdapter is friend of PreciceBC
class PreciceAdapter;

class BoundaryCondition
{
    private:
    protected:
        const int bdr_attr;
        mfem::Coefficient* coeff;
        
    public:
        BoundaryCondition(const int attr) : bdr_attr(attr) {};
        int GetBdrAttr() const { return bdr_attr; }
        mfem::Coefficient& GetCoeffRef() const { return *coeff; }; // To be used only for assigning to linear/bilinear forms or projecting coeff, so declared const
                                                                   // ProjectBdrCoefficient requires non-constant Coefficient&, for example
        virtual void UpdateCoeff(const double time) = 0;
        virtual bool IsEssential() const = 0;
        virtual bool IsConstant() const = 0; // true if d/dt is 0 for this coefficient
        virtual std::string GetInitString() const = 0;

        virtual ~BoundaryCondition() {delete coeff;};
};

class UniformConstantBC : public BoundaryCondition
{   
    private:
    protected:
        const double uniform_value;
    public:
        UniformConstantBC(const int attr, const double in_value);
        bool IsConstant() const { return true; };
        double GetValue() const { return uniform_value; };
        void UpdateCoeff(const double time) {};
        virtual std::string GetInitString() const = 0;
};

class UniformConstantDirichletBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformConstantDirichletBC(const int attr, const double const_value) : UniformConstantBC(attr, const_value){};
        bool IsEssential() const { return true; };
        std::string GetInitString() const;
};

class UniformConstantNeumannBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformConstantNeumannBC(const int attr, const double const_value) : UniformConstantBC(attr, const_value){};
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

    public:
        UniformSinusoidalBC(const int attr, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert);
        bool IsConstant() const { return false; }
        void UpdateCoeff(const double time);
        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;

};

class UniformSinusoidalDirichletBC : public UniformSinusoidalBC
{   
    private:
    protected:
    public:
        UniformSinusoidalDirichletBC(const int attr, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert) : UniformSinusoidalBC(attr, in_amp, in_angfreq, in_phase, in_vert) {};
        bool IsEssential() const { return true; }
        std::string GetInitString() const;

};

class UniformSinusoidalNeumannBC : public UniformSinusoidalBC
{   
    private:
    protected:
    public:
        UniformSinusoidalNeumannBC(const int attr, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert) : UniformSinusoidalBC(attr, in_amp, in_angfreq, in_phase, in_vert) {};
        bool IsEssential() const { return false; }
        std::string GetInitString() const;

};

class PreciceBC : public BoundaryCondition
{
    friend class JOTSSolverInterface; // Allow JOTSSolverInterface access to everything

    private:
    protected:
        const std::string mesh_name;
        const double default_value;
        const std::string read_data_name;
        const std::string write_data_name;
        const int dim;
        bool update_flag;

        int num_dofs;

        double* coords;
        double* read_data_arr;
        double* write_data_arr;

        // Set by adapter:
        int mesh_id;
        int* vertex_ids;
        int read_data_id;
        int write_data_id;

        mfem::Array<int> bdr_elem_indices;
        mfem::Array<int> bdr_dof_indices; // Includes non-true DOFs
        mfem::ParGridFunction coeff_gf;

    public:
        PreciceBC(const int attr, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const double in_value, const std::string in_read, const std::string in_write);
        bool IsConstant() const { return false; };
        void UpdateCoeff(const double time);

        virtual void RetrieveWriteData(const mfem::ParGridFunction &u_gf) = 0;


        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
        ~PreciceBC();
       
};


class PreciceIsothermalBC : public PreciceBC
{
    private:
        const MaterialProperty& k_prop; // Retain k model for sending HF
    protected:
    public:
        PreciceIsothermalBC(const int attr, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const double in_value, const MaterialProperty &in_k) : PreciceBC(attr, f, in_mesh, in_value, "Temperature", "Heat-Flux"), k_prop(in_k) {};
        void RetrieveWriteData(const mfem::ParGridFunction &u_gf);
        std::string GetInitString() const;

        bool IsEssential() const { return true; };
};

class PreciceHeatFluxBC : public PreciceBC
{
    private:
    protected:
    public:
        PreciceHeatFluxBC(const int attr, mfem::ParFiniteElementSpace& f, const std::string in_mesh, const double in_value) : PreciceBC(attr, f, in_mesh, in_value, "Heat-Flux", "Temperature") {};
        void RetrieveWriteData(const mfem::ParGridFunction &u_gf);
        std::string GetInitString() const;

        bool IsEssential() const { return false; };
};