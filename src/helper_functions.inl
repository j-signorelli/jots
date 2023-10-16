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