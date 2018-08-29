/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   05.08.2018
 */
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <vector>
#include <complex>
#include <PlotPy.hpp>



#include "Exception.hpp"

/**
 * Load the root finders themselves
 */
#include "RootBisection.hpp"
#include "RootInterpolation.hpp"
#include "RootSecant.hpp"
#include "RootBrent.hpp"
#include "RootNewtonRaphson.hpp"
#include "RootRidder.hpp"

#include "Allocator.hpp"

namespace anpi
{
// namespace encapsulating benchmarking functionality
namespace bm
{

  //plot points per function
  
  struct plotPoints{
    std::size_t funcNUm = 4;
      std::vector<double> epsilons;
       std::vector<double> functionCalls;
  
      std::vector<double> getEpsilons(){return epsilons;};
      std::vector<double> getFunctionCalls(){return functionCalls;};
  };

/// Square of a number
template <typename T>
inline T sqr(const T x) { return x * x; }

/// Cube of a number
template <typename T>
inline T cube(const T x) { return x * x * x; }

/// First testing function for roots |x|=e^(-x)
template <typename T>
T t1(const T x) { return std::abs(x) - std::exp(-x); }

/// Second testing function for roots e^(-x²) = e^(-(x-3)²/3 )
template <typename T>
T t2(const T x) { return std::exp(-x * x) - std::exp(-sqr(x - T(3)) / T(3)); }

/// Third testing function for roots x² = atan(x)
template <typename T>
T t3(const T x) { return x * x - std::atan(x); }

/// Fourth testing function for roots (x-2)⊃3; + 0.01(x-2)
template <typename T>
T t4(const T x)
{
  const T x0 = x - T(2);
  return cube(x0) + T(0.01) * x0;
}

/**
     * Wrapper class to count function calls
     *
     * This wrapper fulfills the requirements to act as a
     * std::function<T(T)>, and it simply counts the number
     * of calls made to the operator(), before calling
     * the functor provided at construction time.
     */
template <typename T>
class CallCounter
{
protected:
  /// Maximum allowed size for the square matrices
  mutable size_t _counter;

  std::function<T(T)> _f;

public:
  /// Construct
  CallCounter(std::function<T(T)> f) : _counter(0u), _f(f) {}

  /// Access the counter
  inline size_t counter() const { return _counter; }

  /// Reset the counter
  inline void reset() { _counter = 0u; }

  /// Call the function
  T operator()(const T x)
  {
    ++_counter;
    return _f(x);
  }

  /// Call the function
  T operator()(const T x) const
  {
    ++_counter;
    return _f(x);
  }
};

/**
     * Test the given _closed_ root finder
     *
     * The solver must be itself a std::function expecting another
     * std::function (the one whose roots are being looked for), the
     * two limits of the interval enclosing the root, and the
     * tolerance.
     *
     * The tolerances will start from "start", then progressing with
     *   eps = eps*factor
     * until the end value is reached.
     */
template <typename T>
void rootBench(const std::function<T(const std::function<T(T)> &,
                                     T,
                                     T,
                                     const T)> &solver,
               const T start,
               const T end, 
               const T factor)
{

  if ((factor >= static_cast<T>(1)) &&
      (factor < static_cast<T>(0)))
  {
    throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
  }

  // Alias of the function type, for which the roots are being looked for.
  typedef std::function<T(T)> f_type;

  // Try a series of tolerances
  for (T eps = start; eps > end; eps *= factor)
  {
    std::cout << "eps=" << eps << "; ";

    // Create an std::function instance, which wraps the function
    // t1 with the function counter
    f_type c1(CallCounter<T>(t1<T>));
    solver(c1, T(0), T(2), eps);
    std::cout << c1.template target<CallCounter<T>>()->counter() << "; ";

    // now the same with function t2
    f_type c2(CallCounter<T>(t2<T>));
    solver(c2, T(0), T(2), eps);
    std::cout << c2.template target<CallCounter<T>>()->counter() << "; ";

    // now the same with function t3
    f_type c3(CallCounter<T>(t3<T>));
    solver(c3, T(0), T(0.5), eps);
    std::cout << c3.template target<CallCounter<T>>()->counter() << "; ";

    // now the same with function t4
    f_type c4(CallCounter<T>(t4<T>));
    solver(c4, T(1), T(3), eps);
    std::cout << c4.template target<CallCounter<T>>()->counter() << std::endl;
  }
}

/**
     * Test the given _open_ root finder
     *
     * The solver must be itself a std::function expecting another
     * std::function (the one whose roots are being looked for), the
     * starting root guess, and the tolerance.
     */
template <typename T>
void rootBench(const std::function<T(const std::function<T(T)> &,
                                     T,
                                     const T)> &solver,
               const T start,
               const T end,
               const T factor)
{

  if ((factor >= static_cast<T>(1)) &&
      (factor < static_cast<T>(0)))
  {
    throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
  }

  // Alias of the function type, for which the roots are being looked for.
  typedef std::function<T(T)> f_type;

  // Try a series of tolerances
  for (T eps = start; eps > end; eps *= factor)
  {
    std::cout << "eps=" << eps << "; ";

    // Create an std::function instance, which wraps the function
    // t1 with the function counter
    f_type c1(CallCounter<T>(t1<T>));
    solver(c1, T(0), eps);
    std::cout << c1.template target<CallCounter<T>>()->counter() << "; ";

    // now the same with function t2
    f_type c2(CallCounter<T>(t2<T>));
    solver(c2, T(2), eps);
    std::cout << c2.template target<CallCounter<T>>()->counter() << "; ";

    // now the same with function t3
    f_type c3(CallCounter<T>(t3<T>));
    solver(c3, T(0), eps);
    std::cout << c3.template target<CallCounter<T>>()->counter() << "; ";

    // now the same with function t4
    f_type c4(CallCounter<T>(t4<T>));
    solver(c4, T(1), eps);
    std::cout << c4.template target<CallCounter<T>>()->counter() << std::endl;
  }
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*****OVERLOADED FUNCTIONS FOR PLOTTING*****/
/**
     * Test the given _closed_ root finder
     *
     * The solver must be itself a std::function expecting another
     * std::function (the one whose roots are being looked for), the
     * two limits of the interval enclosing the root, and the
     * tolerance.
     *
     * The tolerances will start from "start", then progressing with
     *   eps = eps*factor
     * until the end value is reached.
     */
template <typename T>
void rootBench(const std::function<T(const std::function<T(T)> &,
                                     T,
                                     T,
                                     const T)> &solver,
               const T start,
               const T end, 
               const T factor,
               std::vector<plotPoints>& points)
{

  if ((factor >= static_cast<T>(1)) &&
      (factor < static_cast<T>(0)))
  {
    throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
  }

  // Alias of the function type, for which the roots are being looked for.
  typedef std::function<T(T)> f_type;

  //vector with the 4 pairs of plt point per each test function
plotPoints plot_t1, plot_t2, plot_t3, plot_t4;

//std::vector<plotPoints> points= {plot_t1, plot_t2, plot_t3, plot_t4};
//variable to hold the function call count
size_t functCallCount =0;
//counter to add the poitns


  // Try a series of tolerances
  for (T eps = start; eps > end; eps *= factor)
  {
    std::cout << "eps=" << eps << "; ";
    //add the epsilon to our plot points


    // Create an std::function instance, which wraps the function
    // t1 with the function counter
    f_type c1(CallCounter<T>(t1<T>));
    solver(c1, T(0), T(2), eps);
    functCallCount =c1.template target<CallCounter<T>>()->counter();
    std::cout <<  functCallCount << "; ";
    plot_t1.epsilons.push_back(eps);
    plot_t1.functionCalls.push_back(functCallCount);
    


  //   // now the same with function t2
    f_type c2(CallCounter<T>(t2<T>));
    solver(c2, T(0), T(2), eps);
    functCallCount =c2.template target<CallCounter<T>>()->counter();
    std::cout <<  functCallCount << "; ";
  plot_t2.epsilons.push_back(eps);
    plot_t2.functionCalls.push_back(functCallCount);

  //   // now the same with function t3
  //   f_type c3(CallCounter<T>(t3<T>));
  //   solver(c3, T(0), T(0.5), eps);
  //   functCallCount =c3.template target<CallCounter<T>>()->counter();
  //   std::cout <<  functCallCount << "; ";
   
  //  points[2].getEpsilons().insert(eps);
  //   points[2].getFunctionCalls().insert((T)functCallCount);

  //   // now the same with function t4
  //   f_type c4(CallCounter<T>(t4<T>));
  //   solver(c4, T(1), T(3), eps);
  //   functCallCount =c4.template target<CallCounter<T>>()->counter();
  //   std::cout <<  functCallCount << "; ";
  //   points[3].getEpsilons().insert(eps);
  //   points[3].getFunctionCalls().insert((T)functCallCount);
  }

points.push_back(plot_t1);
points.push_back(plot_t2);
}

/**
     * Test the given _open_ root finder
     *
     * The solver must be itself a std::function expecting another
     * std::function (the one whose roots are being looked for), the
     * starting root guess, and the tolerance.
     */
template <typename T>
void rootBench(const std::function<T(const std::function<T(T)> &,
                                     T,
                                     const T)> &solver,
               const T start,
               const T end,
               const T factor,
               std::vector<plotPoints>& functpoints)
{

  if ((factor >= static_cast<T>(1)) &&
      (factor < static_cast<T>(0)))
  {
    throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
  }

  // Alias of the function type, for which the roots are being looked for.
  typedef std::function<T(T)> f_type;
 //vector with the 4 pairs of plt point per each test function
std::vector<plotPoints> points(4);
//variable to hold the function call count
size_t functCallCount =0;

  // Try a series of tolerances
  for (T eps = start; eps > end; eps *= factor)
  {
    std::cout << "eps=" << eps << "; ";

    // Create an std::function instance, which wraps the function
    // t1 with the function counter
  //   f_type c1(CallCounter<T>(t1<T>));
  //   solver(c1, T(0), eps);
  //    functCallCount =c1.template target<CallCounter<T>>()->counter();
  //   std::cout <<  functCallCount << "; ";
  //   points[0].getEpsilons().insert(eps);
  //   points[0].getFunctionCalls().insert((T)functCallCount);
    

  //   // now the same with function t2
  //   f_type c2(CallCounter<T>(t2<T>));
  //   solver(c2, T(2), eps);
  //    functCallCount =c2.template target<CallCounter<T>>()->counter();
  //   std::cout <<  functCallCount << "; ";
  //   points[1].getEpsilons().insert(eps);
  //   points[1].getFunctionCalls().insert((T)functCallCount);

  //   // now the same with function t3
  //   f_type c3(CallCounter<T>(t3<T>));
  //   solver(c3, T(0), eps);
  //    functCallCount =c3.template target<CallCounter<T>>()->counter();
  //   std::cout <<  functCallCount << "; ";
  //  points[2].getEpsilons().insert(eps);
  //   points[2].getFunctionCalls().insert((T)functCallCount);

  //   // now the same with function t4
  //   f_type c4(CallCounter<T>(t4<T>));
  //   solver(c4, T(1), eps);
  //   functCallCount =c4.template target<CallCounter<T>>()->counter();
  //   std::cout <<  functCallCount << "; ";
  //   points[3].getEpsilons().insert(eps);
  //   points[3].getFunctionCalls().insert((T)functCallCount);
  }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////










/**
     * Benchmark all solvers using a range of tolerances geometrically changing
     * multiplying from the start point until the end with the given factor
     */
template <typename T>
void allSolvers(const T start, const T end, const T factor)
{

  std::cout << "Bisection" << std::endl;
  anpi::bm::rootBench<T>(anpi::rootBisection<T>, start, end, factor);

  std::cout << "Interpolation" << std::endl;
  anpi::bm::rootBench<T>(anpi::rootInterpolation<T>, start, end, factor);

  std::cout << "Secant" << std::endl;
  anpi::bm::rootBench<T>(anpi::rootSecant<T>, start, end, factor);

  std::cout << "NewtonRaphson" << std::endl;
  anpi::bm::rootBench<T>(anpi::rootNewtonRaphson<T>, start, end, factor);

  std::cout << "Brent" << std::endl;
  anpi::bm::rootBench<T>(anpi::rootBrent<T>, start, end, factor);

  std::cout << "Ridder" << std::endl;
  anpi::bm::rootBench<T>(anpi::rootRidder<T>, start, end, factor);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
     * Benchmark all solvers using a range of tolerances geometrically changing
     * multiplying from the start point until the end with the given factor
     */
template <typename T>
void allSolversPlotted(const T start, const T end, const T factor)
{
  //initialize plotting
static anpi::Plot2d<double> plotT1, plotT2, plotT3, plotT4;
plotT1.initialize(1);
plotT1.setTitle("Funcion T1");
plotT2.initialize(2);
plotT2.setTitle("Funcion T2");


std::vector<plotPoints> functPoints;

  std::cout << "Bisection" << std::endl;
  anpi::bm::rootBench<T>(anpi::rootBisection<T>, start, end, factor, functPoints);


  std::cout << "the value of the first epsilon for plotting is:" << functPoints[0].epsilons[0]<< std::endl;
  std::cout << "the value of the first call functions count for plotting is:" << functPoints[0].functionCalls[0]<< std::endl;

std::cout << std::endl;
std::cout <<"tamano de epsilons"<< functPoints[0].epsilons.size();
std::cout << std::endl;


  //for test function t1
 plotT1.plot(functPoints[0].epsilons,functPoints[0].functionCalls,"Biseccion","red");
plotT2.initialize(2);
 plotT2.plot(functPoints[1].epsilons,functPoints[1].functionCalls,"Biseccion","red");

 plotT1.show();
 plotT2.show();

  // std::cout << "Interpolation" << std::endl;
  // anpi::bm::rootBench<T>(anpi::rootInterpolation<T>, start, end, factor, functPoints);

  // std::cout << "Secant" << std::endl;
  // anpi::bm::rootBench<T>(anpi::rootSecant<T>, start, end, factor, functPoints);

  // std::cout << "NewtonRaphson" << std::endl;
  // anpi::bm::rootBench<T>(anpi::rootNewtonRaphson<T>, start, end, factor, functPoints);

  // std::cout << "Brent" << std::endl;
  // anpi::bm::rootBench<T>(anpi::rootBrent<T>, start, end, factor, functPoints);

  // std::cout << "Ridder" << std::endl;
  // anpi::bm::rootBench<T>(anpi::rootRidder<T>, start, end, factor, functPoints);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace bm
} // namespace anpi

BOOST_AUTO_TEST_SUITE(RootFinders)


/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE(RootFinders)
{

  // Benchmark the solvers using float
  std::cout << "<float>" << std::endl;
  anpi::bm::allSolvers<float>(0.1f, 1.e-7f, 0.125f);

  // Benchmark the solvers using double
  std::cout << "<double>" << std::endl;
  anpi::bm::allSolvers<double>(0.1f, 1.e-15f, 0.125f);
}



BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(RootFindersPlotted)

BOOST_AUTO_TEST_CASE(RootFindersPlotted)
{

  // Benchmark the solvers using float
  std::cout << "<float>" << std::endl;
  anpi::bm::allSolversPlotted<float>(0.1f, 1.e-7f, 0.125f);

  // // Benchmark the solvers using double
  // std::cout << "<double>" << std::endl;
  // anpi::bm::allSolversPlotted<double>(0.1f, 1.e-15f, 0.125f);
}

BOOST_AUTO_TEST_SUITE_END()