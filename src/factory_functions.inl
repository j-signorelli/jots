using namespace mfem;

namespace Factory
{

inline IterativeSolver* GetSolver(std::string solver_label, MPI_Comm comm_)
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

inline HypreSmoother::Type GetPrec(std::string prec_label) 
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