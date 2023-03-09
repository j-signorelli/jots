#include <string>
#include <iostream>
using namespace std;

class Config
{
    private:
        string m_input_file;

        string m_mesh_file;
        int m_fe_order;
        int m_serial_refine;
        int m_parallel_refine;
        
        double m_kappa;

        int m_time_scheme;
        double m_tf;
        double m_dt;

        int m_restart_freq;
        int m_vis_freq;
        
    public:
        Config(const char* input_file);

        string GetMeshFile() const 
        {
            return m_mesh_file;
        }

        int GetFEOrder() const
        {
            return m_fe_order;
        }

        int GetTimeIntegration() const
        {
            return m_time_scheme;
        }
        
        double GetFinalTime() const
        {
            return m_tf;
        }

        double Getdt() const
        {
            return m_dt;
        }

        double GetKappa() const
        {
            return m_kappa;
        }

    protected:

};