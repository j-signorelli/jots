using namespace mfem;
using namespace std;

namespace Factory
{

inline IterativeSolver* GetSolver(string solver_label, MPI_Comm comm_)
{
    switch (Linear_Solver_Map.at(solver_label))
    {
        case LINEAR_SOLVER::CG:
            return new CGSolver(comm_);
            break;
        case LINEAR_SOLVER::GMRES:
            return new GMRESSolver(comm_);
            break;
        case LINEAR_SOLVER::FGMRES:
            return new FGMRESSolver(comm_);
            break;
    }

    return nullptr;
}

inline HypreSmoother::Type GetPrec(string prec_label) 
{
    switch (Preconditioner_Map.at(prec_label))
    {
        case PRECONDITIONER::JACOBI:
            return HypreSmoother::Jacobi;
            break;
        case PRECONDITIONER::CHEBYSHEV:
            return HypreSmoother::Chebyshev;
            break;
    }

    return HypreSmoother::Jacobi;
}

inline ODESolver* GetODESolver(string time_scheme_label)
{
    switch (Time_Scheme_Map.at(time_scheme_label))
    {
        case TIME_SCHEME::EULER_IMPLICIT:
            return new BackwardEulerSolver;
            break;
        case TIME_SCHEME::EULER_EXPLICIT:
            return new ForwardEulerSolver;
            break;
        case TIME_SCHEME::RK4:
            return new RK4Solver;
    }

    return nullptr;
}

inline IterativeSolver::PrintLevel CreatePrintLevel(std::vector<std::string> print_level_settings)
{
    IterativeSolver::PrintLevel pl;
    for (size_t i = 0; i < print_level_settings.size(); i++)
    {
        switch (Print_Level_Map.at(print_level_settings[i]))
        {
            case PRINT_LEVEL::NONE:
                pl.None();
                break;
            case PRINT_LEVEL::WARNINGS:
                pl.Warnings();
                break;
            case PRINT_LEVEL::ERRORS:
                pl.Errors();
                break;
            case PRINT_LEVEL::ITERATIONS:
                pl.Iterations();
                break;
            case PRINT_LEVEL::FIRSTANDLAST:
                pl.FirstAndLast();
                break;
            case PRINT_LEVEL::SUMMARY:
                pl.Summary();
                break;
            case PRINT_LEVEL::ALL:
                pl.All();
                break;
            default:
                break;
        }
    }
    return pl;

}


}

namespace Helper
{
    
template<typename Key, typename Value>
inline vector<Key> GetKeyVector(map<Key, Value> in_map)
{
    vector<Key> keys;
    for (typename map<Key, Value>::iterator it = in_map.begin(); it != in_map.end(); it++)
    {
        keys.push_back(it->first);
    }

    return keys;
}

}

inline void JOTSNewtonSolver::SetOperator(const Operator &op)
{
    NewtonSolver::SetOperator(op);

    // Check if iterating on k
    if (dynamic_cast<const JOTS_k_Operator*>(&op) != nullptr)
        iterate_on_k = true;
    else
        iterate_on_k = false;
}

inline void JOTSNewtonSolver::ProcessNewState(const Vector& x) const
{
    for (int i = 0; i < mps.Size(); i++)
    {
        if (!mps[i]->IsConstant())
        {
            if (iterate_on_k) // x = k
            {   
                // Retrieve u_n and dt from the operator
                const JOTS_k_Operator* jop = dynamic_cast<const JOTS_k_Operator*>(oper);
                const Vector& u_n = jop->Get_u_n();
                const double& dt = jop->Get_dt();
                
                // Update coefficients
                Vector z(x.Size());
                add(u_n, dt, x, z);
                mps[i]->UpdateAllCoeffs(z);
            
            }
            else // x = u
                mps[i]->UpdateAllCoeffs(x);
        }
    }
}