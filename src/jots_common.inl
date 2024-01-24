using namespace mfem;
using namespace std;

namespace Factory
{

inline IterativeSolver* GetSolver(string solver_label, MPI_Comm comm_)
{
    switch (Solver_Map.at(solver_label))
    {
        case SOLVER::CG:
            return new CGSolver(comm_);
            break;
        case SOLVER::GMRES:
            return new GMRESSolver(comm_);
            break;
        case SOLVER::FGMRES:
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

void JOTSNewtonSolver::SetOperator(const Operator& op) override
{
    NewtonSolver::SetOperator(op);
    // If this is a JOTS_k_Operator, update flag.
    if (dynamic_cast<JOTS_k_Operator*>(&op) != nullptr)
    {
        iterate_on_k = true;
    }
    else
        iterate_on_k = false;
}

void JOTSNewtonSolver::ProcessNewState(const Vector& x) override
{
    for (int i = 0; i < mps.Size(); i++)
    {
        if (!mps[i]->IsConstant())
        {
            if (iterate_on_k) // x = k
            {   
                // Retrieve u_n and dt from the operator
                const Vector& u_n = dynamic_cast<JOTS_k_Operator>(op).Get_u_n();
                const double& dt = dynamic_cast<JOTS_k_Operator>(op).Get_dt();
                
                // Update coefficients
                mps[i].UpdateAllCoeffs(u_n + dt*x);
            
            }
            else // x = u
                mp[i].UpdateAllCoeffs(x);
        }
    }
}