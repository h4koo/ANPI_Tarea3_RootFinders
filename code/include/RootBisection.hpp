/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the bisection method.
   *
   * @param funct a std::function of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps) {
    
    //not sure if we want this here
const int MAX_ITERATIONS=40;
    // TODO: Put your code in here!
//evaluate lower boundary
  T f_min = funct(xl);
  int iterations=0;
  //while the difference between the boundaries is more than the desired accuracy
    while (xl + eps < xu) {

      //calculate mid value and evaluate function
        T const mid = 0.5 * xl + 0.5 * xu;
        T const f_mid = funct(mid);

        //if both our mid value and min value are negative
        //then our mid value is our new lower
        if ((f_min < 0) && (f_mid < 0)) {
            xl = mid;
            f_min = f_mid;
        } else {//our new upper would me the mid value
            xu = mid;
        }
        ++iterations;
        if (iterations==MAX_ITERATIONS) return std::numeric_limits<T>::quiet_NaN();
    }

    return xl;
    // Return NaN if no root was found 
}//end if sign of boundaries is different





  }


  
#endif

