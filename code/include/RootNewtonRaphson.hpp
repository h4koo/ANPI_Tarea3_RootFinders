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

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking by means of the
   * Newton-Raphson method
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial root guess
   * 
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>   
  T dfunct(const std::function<T(T)>& funct,T xi){
	  double h = 0.00001;
	  double derivada = 0;
	  derivada= (funct(xi+h)-funct(xi))/h;
	  return derivada;

  }

  
  template<typename T>
  T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps) {

    // TODO: Put your code in here!
    T xii;
	  T Dx;
	
	  do {
		  xii = xi - funct(xi)/dfunct(funct,xi);
		  Dx = fabs(xii - xi);
		  xi = xii;
	  } while (Dx > eps);
	  return xii;

    // Return NaN if no root was found
     //return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif
