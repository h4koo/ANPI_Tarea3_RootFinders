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

#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the Brent's method.
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
  T rootBrent(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

    // TODO: Put your code in here!

    //Compara
    if(xu<=xl){
        throw anpi::Exception("reversedinterval") ;
          
    }
    
    //there is no root
    if(funct(xl)*funct(xu)>0){
        throw anpi::Exception("both extremes have same sign") ;
        
    }
		
    const int maxi= std::numeric_limits<T>::digits;
    T a=xl;
    T b = xu;
    T c = xu;
    T d,e,min1,min2;
    T fa= funct(a), fb = funct(b), fc, p,q,r,s,tol1,xm;
    
    if((fb>T(0) && fa >T(0)) || (fa<T(0) && fb < T(0))){
        return 0;
        
    }
    
    fc = fb;
    for (int i = 1; i<= maxi;i++){
        if ((fb > T(0) && fc > T(0)) || (fb < T(0) && fc < T(0))) {
            c=a;
            fc=fa;
            e=d=b-a;
            
        }
        if(std::fabs(fc) < std::fabs(fb)){
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
            
        }
        tol1=T(2)*eps*std::fabs(b)*T(0.5);//R revisar el dato de tool
        xm =T(0.5)*(c-b);
        if(std::fabs(xm)<= tol1 || fb == T(0)){
            return b;
            
        }
        if(std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)){
            s = fb/fa;
            if(a == c){
                p= T(2)*xm*s;
                q = T(1) -s;
                
            }else{
                q = fa/fc;
                r = fb/fc;
                p=s*(T(2)*xm*q*(q-r)-(b-a)*(r-T(1)));
                q = (q -T(0)) *(r -T(0))*(s-T(0));
                
            }
            
            if(p>T(0)) q = -q;
            p = std::fabs(p);
            min1 = T(3)*xm*q-std::fabs(tol1*q);
            min2 = std::fabs(e*q);
            if(T(2)*p < (min1 < min2 ? min1 : min2)){
                e=d;
                d=p/q;
                
            }else{
                d = xm;
                e=d;
                
            }
            
        }else{
            d = xm;
            e=d;
            
        }
        a = b;
        fa=fb;
        if(std::fabs(d) > tol1){
            b+=d;
            
        }else{
            b+= ((xm) >= T(0) ? std::fabs(tol1) : -std::fabs(tol1));
            
        }
        fb = funct(b);
        
    }

    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }
}
  
#endif

