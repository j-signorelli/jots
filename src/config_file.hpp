#include <string>

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

        const char* GetMeshFile()
        {
            return mesh_file;
        }

        const int GetFEOrder()
        {
            return fe_order;
        }

        const int GetTimeIntegration()
        {
            return time_integration;
        }
        
        const double GetFinalTime()
        {
            return tf;
        }

        const double Getdt()
        {
            return dt;
        }

    protected:

};