#pragma once

class ThermDiff
{
    private:
    
    public:

    protected:
};

class ConstantThermDiff : public ThermDiff
{
    private:
        double kappa;

    public:
        ConstantThermDiff(double in_kappa);

    protected:
};