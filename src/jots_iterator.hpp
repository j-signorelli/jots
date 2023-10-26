#pragma once

#include "mfem/mfem.hpp"

#include "option_structure.hpp"

class JOTSIterator
{
    private:
    protected:
        static constexpr double TIME_TOLERANCE = 1e-14;
    public:
        virtual bool IsNotComplete() const = 0;
        virtual void Iterate(mfem::Vector& u) = 0;
        virtual void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) = 0;
        virtual void UpdateNeumann() = 0;
        virtual ~JOTSIterator() = default;
};