#include <string>

/*!
 * This idea below was taken from SU2 v6.0:
 *
 * \class CCreateMap
 * \brief creates a map from a list by overloading operator()
 * \tparam T - the key type in the map
 * \tparam U - the mapped value type in the map
 * \author Boost.Assign and anonymous person on stackoverflow
 *
 * We need this to create static const maps that map strings to enum
 * types.  The implementation is based on the Boost.Assign library.  This
 * particular version is taken from
 * http://stackoverflow.com/questions/138600/initializing-a-static-stdmapint-int-in-c
 */
template <typename T, typename U>
class CCreateMap
{
 private:
  std::map<T, U> m_map;

 public:
  CCreateMap(const T& key, const U& val)
  {
    m_map[key] = val;
  }
  CCreateMap<T, U>& operator()(const T& key, const U& val)
  {
    m_map[key] = val;
    return *this;
  }
  operator std::map<T, U>()
  {
    return m_map;
  }
};

enum class TIME_INTEGRATION
{

    // Implicit L-stable methods
    BACKWARD_EULER = 0,
    //case 2:  ode_solver = new SDIRK23Solver(2); break;
    //case 3:  ode_solver = new SDIRK33Solver; break;

    // Explicit methods
    FORWARD_EULER = 1,
    RK4 = 2
    //case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
    //case 13: ode_solver = new RK3SSPSolver; break;
    //case 14: ode_solver = new RK4Solver; break;
    //case 15: ode_solver = new GeneralizedAlphaSolver(0.5); break;

    // Implicit A-stable methods (not L-stable)
    //case 22: ode_solver = new ImplicitMidpointSolver; break;
    //case 23: ode_solver = new SDIRK23Solver; break;
    //case 24: ode_solver = new SDIRK34Solver; break;
}

static const map<string, TIME_INTEGRATION> Time_Integration_Map = CCreateMap<string, TIME_INTEGRATION>("BACKWARD_EULER", TIME_INTEGRATION::BACKWARD_EULER)("FORWARD_EULER", TIME_INTEGRATION::FORWARD_EULER)("RK4", TIME_INTEGRATION::RK4)