#include "option_structure.hpp"

class Config
{
    private:
        string input_file;
        char* mesh_file;
        int fe_order;
        int time_integration;
        double tf;
        double dt;
        double kappa;
        double alpha;
        int output_restart_freq;
        int output_vis_freq;
        
    public:
        Config(string in_file);

        char* GetMeshFile()
        {
            return mesh_file;
        }

        int GetFEOrder()
        {
            return fe_order;
        }

        int GetTimeIntegration()
        {
            return time_integration;
        }
        
        double GetFinalTime()
        {
            return tf;
        }

        double Getdt()
        {
            return dt;
        }

    protected:

};