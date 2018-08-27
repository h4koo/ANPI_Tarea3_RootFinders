/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 04.08.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_RIDDER_HPP
#define ANPI_ROOT_RIDDER_HPP

namespace anpi
{

/**
   * Find a root of the function funct looking for it starting at xi
   * by means of the secant method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial position
   * @param xii second initial position 
   *
   * @return root found, or NaN if no root could be found
   */
template <typename T>
T rootRidder(const std::function<T(T)> &funct, T xi, T xii, const T eps)
{

  // TODO: Put your code in here!
  //max amount of iterations before the function returns NaN
  const int MAX_ITERATION = 40;

  //values of function we want to look the roots for at the provided
  //boundaries
  T fl = funct(xi);
  T fh = funct(xii);

  //the value has to be enclosed, therefore, the function evaluations must have diferent sign
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0))
  {
    T xl = xi;
    T xh = xii;
    
    //Initialize with any highly unlikely value, to simplify logic below.**
    T ans = -9.99e99; 

    //variables that will hold values for the calculation in the for loop
    T xm, fm, s, xnew, fnew;

    for (int j = 0; j < MAX_ITERATION; j++)
    {

      //middle point between our current boundaries
      xm = 0.5 * (xl + xh);
      //evaluate function at middle point
      fm = func(xm);

      //s is used to calculate new boundary xnew
      s = sqrt(fm * fm - fl * fh);

      //this means fm==fl==fh, we have reached a solution
      if (s == 0.0)
        return ans;

      //Calculating xnew as per Ridder
      xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);

      //if the difference between our new point and our current answer is less than the desired
      //accuracy (eps) then we have reached a solution. ** This is why we need to initialize ans
      if (abs(xnew - ans) <= eps)
        return ans;

      //set the answer to our estimation
      ans = xnew;
      //evaluate function at estimation
      fnew = func(ans);

      //if the evaluated function at our estimation is 0
      //we have reached a solution
      if (fnew == 0.0)
        return ans;

      //if the middle value and estimated value have different signs
      if ((fm > 0 && fnew < 0) || (fm < 0 && fnew > 0))
      {

        //set our lower boundary to our middle value
        xl = xm;
        fl = fm;

        //set higher boundary to our estimated value
        xh = ans;
        fh = fnew;
      }
      else if ((fl > 0 && fnew < 0) || (fl < 0 && fnew > 0))
      {
        //if the change in sign is between our lower boundary and our estimation the root is between this range
        xh = ans;
        fh = fnew;
      }
      //if the change in sign is between our higher boundary and our estimation the root is between this range
      else if ((fh > 0 && fnew < 0) || (fh < 0 && fnew > 0))
      {

        xl = ans;
        fl = fnew;
      }
      else
        throw("rootRidder should never get here.");

      //if the difference between our boundaries is less than the desired accuracy we reached a solution
      if (abs(xh - xl) <= eps)
        return ans;
    } //end for

    //solution not found in MAX_ITERATION amout of tries
    throw("rootRidder exceeded maximum iterations");
  }
  else
  {
    //check if one of the boundaries is a solution
    if (fl == 0.0)
      return xi;
    if (fh == 0.0)
      return xii;

    //error reached
    throw("root must be bracketed in rootRidder.");
  }

  // Return NaN if no root was found
  return std::numeric_limits<T>::quiet_NaN();
}

} // namespace anpi

#endif
