/*
This implementation of the Nelder-Mead optimization method is based on a code by
Piotr Różański that can be found here: https://github.com/develancer/nelder-mead
I modified it so that it behaves exactly like the implementation found in the R
package dfoptim by Ravi Varadhan.

This code is distributed under the GNU GPL license.

Author: Jasmine Gamblin
*/

#ifndef PTR_NELDER_MEAD_H
#define PTR_NELDER_MEAD_H

#include <array>
#include <climits>
#include <functional>
#include <cstring>
#include <cmath>
#include <algorithm>


template<typename real, int n>
struct nelder_mead_result {
    std::array<real,n> xmin;
    real ymin;
    int icount;
    std::string message;
};


template<typename real, int n>
std::array<real, n> g(std::array<real, n> &x, std::array<real, n> &lower, std::array<real, n> &upper)
{
    std::array<real,n> res;
    for (int j=0; j<n; j++)
    {
        res[j] = atanh(2.*(x[j]-lower[j])/(upper[j]-lower[j])-1.);
    }
    return res;
}

template<typename real, int n>
std::array<real, n> ginv(std::array<real, n> &x, std::array<real, n> &lower, std::array<real, n> &upper)
{
    std::array<real,n> inv;
    for (int j=0; j<n; j++)
    {
        inv[j] = lower[j] + (upper[j]-lower[j])/2.*(1.+tanh(x[j]));
    }
    return inv;
}


template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


template<typename real, int n>
nelder_mead_result<real,n> nelder_mead(
        const std::function<real(const std::array<real,n> &, int)> &fn,
        std::array<real,n> &start,
        std::array<real,n> &lower,
        std::array<real,n> &upper,
        real tol,
        int maxfeval,
        int restartsMax
) {
    real scale, alpha1, alpha2, simplexSize;
    real y2star, ystar, ynew, aux;
    bool happy;
    std::array<real,n> p[n+1];
    std::array<real,n> p2[n+1];
    std::array<real,n> v[n];
    std::array<real,n> x0, pstar, p2star, pbar, pnew, sgrad;
    real y[n+1];
    real y2[n+1];
    real delf[n];
    real diam[n];

    nelder_mead_result<real,n> result;
    result.icount = 0;
    
    // Simplex initialisation
    x0 = g<double,n>(start, lower, upper);
    aux = 0.;
    for (int i=0; i<n; i++) {aux += std::pow(x0[i],2);}
    scale = max(1.0,sqrt(aux));
    alpha1 = scale/(n*sqrt(2.)) * (sqrt(n+1.) + n-1.);
    alpha2 = scale/(n*sqrt(2.)) * (sqrt(n+1.) - 1.);
    
    p[0] = x0;
    y[0] = fn(ginv<double,n>(x0, lower, upper), result.icount+1);
    result.icount++;
    for (int i = 1; i<n+1; i++)
    {
        for (int j = 0; j<n; j++)
        {
            p[i][j] = x0[j]+alpha2;
        }
        p[i][i-1] = x0[i-1]+alpha1;   
        y[i] = fn(ginv<double,n>(p[i], lower, upper), result.icount+1);
        result.icount++;
    }

    // Sort vertices by Y value
    vector<int> seq(n+1);   
    for (int i = 0 ; i < n+1 ; i++) {seq[i] = i;}
    vector<int> order = seq;
    sort(order.begin(), order.end(), [&](const int& a, const int& b) {return (y[a] < y[b]);});

    real rho = 1.;
    real gamma = 0.5;
    real chi = 2.;
    real sigma = 0.5;
    int restarts = 0;
    bool orth = false;

    // initialisation of some variables
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++) {v[i][j] = p[order[i+1]][j]-p[order[0]][j];}
        delf[i] = y[order[i+1]]-y[order[0]];
    }
    real aux2 = 0.; real aux3 = 0.;
    real aux4;
    for (int i=0; i<n; i++)
    {
        aux = 0.; aux4 = 0.;
        for (int j=0; j<n; j++)
        {
            aux += std::pow(v[i][j],2);
            aux2 += abs(v[i][j]);
            aux4 += v[i][j]*delf[j];
        }
        diam[i] = sqrt(aux);
        aux3 += abs(p[order[0]][i]);
        sgrad[i] = aux4;
    }
    simplexSize = aux2/max(1.0, aux3);
    aux = 0.;
    for (int i=0; i<n; i++) {aux += std::pow(sgrad[i],2);}
    real alpha3 = 1e-04*(*std::max_element(diam, diam+n))/sqrt(aux);
    
    
    // Inner loop.
    while (result.icount < maxfeval && y[order[n]]-y[order[0]] > tol
    && restarts < restartsMax && simplexSize > 1e-6) {
        happy = false;

        // Calculate PBAR, the centroid of the simplex vertices
        // excepting the vertex with Y value YHI.
        for (int i = 0; i < n; i++) {
            aux = 0.;
            for (int j = 0; j < n + 1; j++) {
                aux += p[j][i];
            }
            aux -= p[order[n]][i];
            pbar[i] = aux / n;
        }
        // Reflection through the centroid.
        for (int i = 0; i < n; i++) {
            pstar[i] = (1.+rho)*pbar[i] - rho*p[order[n]][i];
        }
        ystar = fn(ginv<double,n>(pstar, lower, upper), result.icount+1);
        result.icount++;

        // b)
        if (y[order[0]] <= ystar && ystar < y[order[n-1]]) {
            happy = true;
            pnew = pstar;
            ynew = ystar;
        }
        // c)
        else if (ystar < y[order[0]]) {
            for (int i = 0; i < n; i++) {
                p2star[i] = (1.+rho*chi)*pbar[i] - rho*chi*p[order[n]][i];
            }
            y2star = fn(ginv<double,n>(p2star, lower, upper), result.icount+1);
            result.icount++;

            // Check extension.
            if (y2star < ystar) {
                pnew = p2star;
                ynew = y2star;
                happy = true;
            } else {
                // Retain reflection
                pnew = pstar;
                ynew = ystar;
                happy = true;
            }
        }
        // d)
        else if (y[order[n-1]] < ystar && ystar < y[order[n]]) {
            for (int i = 0; i < n; i++) {
                p2star[i] = (1.+rho*gamma)*pbar[i] - rho*gamma*p[order[n]][i];
            }
            y2star = fn(ginv<double,n>(p2star, lower, upper), result.icount+1);
            result.icount++;
            
            if (y2star < ystar) {
                pnew = p2star;
                ynew = y2star;
                happy = true;
            }
        }
        // e)
        else {
            for (int i = 0; i < n; i++) {
                p2star[i] = (1.-gamma)*pbar[i] + gamma*p[order[n]][i];
            }
            y2star = fn(ginv<double,n>(p2star, lower, upper), result.icount+1);
            result.icount++;  

            if (y2star < y[order[n]]) {
                pnew = p2star;
                ynew = y2star;
                happy = true;
            }
        }
        // f)
        if (happy) {
            real armtst = 0;
            for (int i=0; i<n; i++)
            {
                armtst += std::pow(sgrad[i],2);
            }
            armtst *= alpha3;
            if ((y[order[n]]-ynew)/(n+1.) < armtst/n) {
                restarts++;
                orth = true;
                real diams = *std::min_element(diam, diam+n);
                happy = false;
                p2[0] = p[order[0]];
                y2[0] = y[order[0]];
                for (int i=0; i<n; i++)
                {
                    for (int j=0; j<n; j++)
                    {
                        p2[i+1][j] = p[order[0]][j];
                    }
                    p2[i+1][i] -= diams*sign(sgrad[i]);
                    y2[i+1] = y[order[i+1]];
                }
                std::memcpy(p, p2, sizeof(p2));
                std::memcpy(y, y2, sizeof(y2));
                order = seq;
            }
            p[order[n]] = pnew;
            y[order[n]] = ynew;
        }
        else if (!happy && restarts < restartsMax)
        {
            if (!orth) 
            {
                orth = true;
            } 
            p2[0] = p[order[0]];
            for (int i=0; i<n; i++)
            {
                for (int j=0; j<n; j++)
                {
                    p2[i+1][j] = p[order[0]][j] - sigma*(p[order[i+1]][j]-p[order[0]][j]);
                }
            }
            std::memcpy(p, p2, sizeof(p2));
            y[0] = y[order[0]];
            for (int i=1; i<n+1; i++) {y[i] = fn(ginv<double,n>(p[i], lower, upper), result.icount+1);}
        }

        // sort vertices
        sort(order.begin(), order.end(), [&](const int& a, const int& b) {return (y[a] < y[b]);});

        // update some variables
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<n; j++) {v[i][j] = p[order[i+1]][j]-p[order[0]][j];}
            delf[i] = y[order[i+1]]-y[order[0]];
        }
        real aux2 = 0.; real aux3 = 0.;
        real aux4;
        for (int i=0; i<n; i++)
        {
            aux = 0.; aux4 = 0.;
            for (int j=0; j<n; j++)
            {
                aux += std::pow(v[i][j],2);
                aux2 += abs(v[i][j]);
                aux4 += v[i][j]*delf[j];
            }
            diam[i] = sqrt(aux);
            aux3 += abs(p[order[0]][i]);
            sgrad[i] = aux4;
        }
        simplexSize = aux2/max(1.0, aux3);
    }
    
    result.xmin = ginv<double,n>(p[order[0]], lower, upper);
    result.ymin = y[order[0]];
    
    if (y[order[n]]-y[order[0]] <= tol || simplexSize <= 1e-06) {
        result.message = "Successful convergence";
    }
    else if (result.icount >= maxfeval) {
        result.message = "Maximum number of fevals exceeded";
    }
    else if (restarts >= restartsMax) {
        result.message = "Stagnation in Nelder-Mead";
    }
    return result;
}

#endif // PTR_NELDER_MEAD_H
