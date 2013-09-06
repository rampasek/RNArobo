/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : Welch Two Sample t-test
 *                (performs the upper, one-sided test)
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>, Martin Kralik <majak47@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

namespace WT {

double criticalValuesAt80[102] = { -1,
1.376,1.061,0.978,0.941,0.92,0.906,0.896,0.889,0.883,0.879,
0.876,0.873,0.87,0.868,0.866,0.865,0.863,0.862,0.861,0.86,
0.859,0.858,0.858,0.857,0.856,0.856,0.855,0.855,0.854,0.854,
0.853,0.853,0.853,0.852,0.852,0.852,0.851,0.851,0.851,0.851,
0.85,0.85,0.85,0.85,0.85,0.85,0.849,0.849,0.849,0.849,
0.849,0.849,0.849,0.849,0.849,0.848,0.848,0.848,0.848,0.848,
0.848,0.848,0.848,0.847,0.847,0.847,0.847,0.847,0.847,0.847,
0.847,0.847,0.847,0.846,0.846,0.846,0.846,0.846,0.846,0.846,
0.846,0.846,0.846,0.846,0.846,0.846,0.846,0.846,0.846,0.846,
0.846,0.846,0.845,0.845,0.845,0.845,0.845,0.845,0.845,0.845,0.842
};
double criticalValuesAt90[102] = { -1,
3.078,1.886,1.638,1.533,1.476,1.44,1.415,1.397,1.383,1.372,
1.363,1.356,1.35,1.345,1.341,1.337,1.333,1.33,1.328,1.325,
1.323,1.321,1.319,1.318,1.316,1.315,1.314,1.313,1.311,1.31,
1.309,1.309,1.308,1.307,1.306,1.306,1.305,1.304,1.304,1.303,
1.303,1.302,1.302,1.301,1.301,1.3,1.3,1.299,1.299,1.299,
1.298,1.298,1.298,1.297,1.297,1.297,1.297,1.296,1.296,1.296,
1.296,1.295,1.295,1.295,1.295,1.295,1.294,1.294,1.294,1.294,
1.294,1.293,1.293,1.293,1.293,1.293,1.293,1.292,1.292,1.292,
1.292,1.292,1.292,1.292,1.292,1.291,1.291,1.291,1.291,1.291,
1.291,1.291,1.291,1.291,1.291,1.29,1.29,1.29,1.29,1.29,1.282
};

double criticalValuesAt95[102] = { -1,
6.314,2.92,2.353,2.132,2.015,1.943,1.895,1.86,1.833,1.812,
1.796,1.782,1.771,1.761,1.753,1.746,1.74,1.734,1.729,1.725,
1.721,1.717,1.714,1.711,1.708,1.706,1.703,1.701,1.699,1.697,
1.696,1.694,1.692,1.691,1.69,1.688,1.687,1.686,1.685,1.684,
1.683,1.682,1.681,1.68,1.679,1.679,1.678,1.677,1.677,1.676,
1.675,1.675,1.674,1.674,1.673,1.673,1.672,1.672,1.671,1.671,
1.67,1.67,1.669,1.669,1.669,1.668,1.668,1.668,1.667,1.667,
1.667,1.666,1.666,1.666,1.665,1.665,1.665,1.665,1.664,1.664,
1.664,1.664,1.663,1.663,1.663,1.663,1.663,1.662,1.662,1.662,
1.662,1.662,1.661,1.661,1.661,1.661,1.661,1.661,1.66,1.66,1.645 
};

double criticalValuesAt975[102] = { -1,
12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,
2.201,2.179,2.16,2.145,2.131,2.12,2.11,2.101,2.093,2.086,
2.08,2.074,2.069,2.064,2.06,2.056,2.052,2.048,2.045,2.042,
2.04,2.037,2.035,2.032,2.03,2.028,2.026,2.024,2.023,2.021,
2.02,2.018,2.017,2.015,2.014,2.013,2.012,2.011,2.01,2.009,
2.008,2.007,2.006,2.005,2.004,2.003,2.002,2.002,2.001,2,
2,1.999,1.998,1.998,1.997,1.997,1.996,1.995,1.995,1.994,
1.994,1.993,1.993,1.993,1.992,1.992,1.991,1.991,1.99,1.99,
1.99,1.989,1.989,1.989,1.988,1.988,1.988,1.987,1.987,1.987,
1.986,1.986,1.986,1.986,1.985,1.985,1.985,1.984,1.984,1.984,1.96
};

double criticalValuesAt99[102] = { -1,
31.821,6.965,4.541,3.747,3.365,3.143,2.998,2.896,2.821,2.764,
2.718,2.681,2.65,2.624,2.602,2.583,2.567,2.552,2.539,2.528,
2.518,2.508,2.5,2.492,2.485,2.479,2.473,2.467,2.462,2.457,
2.453,2.449,2.445,2.441,2.438,2.434,2.431,2.429,2.426,2.423,
2.421,2.418,2.416,2.414,2.412,2.41,2.408,2.407,2.405,2.403,
2.402,2.4,2.399,2.397,2.396,2.395,2.394,2.392,2.391,2.39,
2.389,2.388,2.387,2.386,2.385,2.384,2.383,2.382,2.382,2.381,
2.38,2.379,2.379,2.378,2.377,2.376,2.376,2.375,2.374,2.374,
2.373,2.373,2.372,2.372,2.371,2.37,2.37,2.369,2.369,2.368,
2.368,2.368,2.367,2.367,2.366,2.366,2.365,2.365,2.365,2.364,2.326
};


/** 
 * Welch Two Sample t-test
 * (performs the upper, two-tailed test)
 * 
 * args:
 * unpaired samples @x and @y with possibly different variance
 * table @criticalValues of critical values for Student's t dist. at some significance level
 * 
 * alternative hypothesis: true difference in means is greater than 0
 * null: real_mean(x) <= real_mean(y)
 * alternative: real_mean(x) > real_mean(y)
 * 
 * returns 1 if the null can be rejected in favor of the alternative hyp.,
 *        -1 if the alternative can be rejected in favor of the null hyp.,
 *         0 otherwise
 */
int welchTest(vector<double> &x, vector<double> &y, double* criticalValues){
    double X1, X2; //sample mean
    double sd1, sd2; //unbiased sample variance (squared deviation)
    unsigned int N1, N2; //sample size
    
    N1 = x.size();
    N2 = y.size();
    if(N1 < 2 || N2 < 2) return 0;   
    
    //calculate sample mean
    double sum1=0, sum2=0;
    for(int i=0;i<N1;++i){
        sum1 += x[i];
    }
    for(int i=0;i<N2;++i){
        sum2 += y[i];
    }
    X1 = sum1/N1;
    X2 = sum2/N2;
    
    //calculate unbiased sample variance (squared deviation)
    sum1=sum2=0;
    for(int i=0;i<N1;++i){
        sum1 += (x[i] - X1)*(x[i] - X1);
    }
    for(int i=0;i<N2;++i){
        sum2 += (y[i] - X2)*(y[i] - X2);
    }
    sd1 = sum1/(N1-1);
    sd2 = sum2/(N2-1);
    
    double g1, g2; //intermediate variables
    g1 = sd1 / N1;
    g2 = sd2 / N2;
    
    double t; // statistic t for Welch's t-test 
    t = (X1 - X2) / sqrt(g1 + g2);
    
    double df; //number of degrees of freedom
    //calculate estimate using the Welch-Satterthwaite equation
    df = (g1 + g2)*(g1 + g2);
    df /= g1*g1/(N1-1) + g2*g2/(N2-1);
       
    //round the number of degrees of freedom (to be used for look up in a table of critical values)
    int roundedDF = (int)round(df);
    if (roundedDF <= 0 ) { roundedDF = 1; }
    if (roundedDF > 100) { roundedDF = 101; } //take "infinite" d.f.
    
    if(t >= criticalValues[roundedDF]){
        return 1;
    } else if(-t >= criticalValues[roundedDF]){
        return -1;
    } else {
        return 0;
    } 
}

/**
 * Make test of Welch's t-test on two examples
 */
bool test(void){
    /*
     * t.test(xx, yy, alternative = "g", paired=FALSE, var.equal=FALSE)
     * 
     *    Welch Two Sample t-test
     * 
     * data:  xx and yy 
     * t = 1.1978, df = 60.947, p-value = 0.1178
     * alternative hypothesis: true difference in means is greater than 0 
     * 95 percent confidence interval:
     * -47.58201       Inf 
     * sample estimates:
     * mean of x mean of y 
     * 2095.676  1975.055 
     * 
     */
    double x1[20] = {
        2224.779,2588.110,1979.625,2137.442,2565.818,1754.023,1654.947,1789.256
        ,2320.659,2039.532,1983.497,2232.903,2513.930,2066.382,2492.715,1988.287
        ,1840.036,2249.749,1766.982,1724.840
    };
    double y1[50] = {
        2465.0984,1909.0328,1175.8747,2171.3780,2193.2821,2854.9475,2060.1777
        ,2258.2366,1856.0535,1501.8126,2987.6542,1681.9778,2479.6776,1259.8584
        ,1120.9043,1982.1213,3012.3949,2252.3730,2591.3122,1940.5890,1995.1850
        ,2535.1344,597.3155,2343.2192,3154.8400,1125.1966,1227.8842,1692.8050
        ,2539.6772,1936.1927,1783.7795,1703.4384,2077.1940,1614.4071,2360.0365
        ,1619.2781,2033.5109,2333.7834,2144.0485,2583.8709,1116.7213,1601.9383
        ,1570.0431,1963.0777,1639.2533,2277.5223,1991.9286,2044.3338,1794.4781,1597.9119
    };
    
    /*
     * t.test(xx, yy, alternative = "g", paired=FALSE, var.equal=FALSE)
     * 
     *    Welch Two Sample t-test
     * 
     * data:  xx and yy 
     * t = 2.0004, df = 79.056, p-value = 0.02445
     * alternative hypothesis: true difference in means is greater than 0 
     * 95 percent confidence interval:
     * 4.351577      Inf 
     * sample estimates:
     * mean of x mean of y 
     * 151.8048  125.8999 
     * 
     */    
    double x2[50] = {
        88.03913,163.31985,162.85867,104.24684,159.79714,219.51456,137.06262
        ,118.78998,132.37481,198.59919,183.85814,114.64763,120.82348,177.12229
        ,96.82845,182.24798,235.88351,86.23593,80.60078,187.91515,119.62568
        ,121.41153,89.96165,165.76889,121.90483,176.79409,94.08487,167.83406
        ,108.77790,240.88811,203.64804,217.18763,158.04343,163.71723,215.27359
        ,121.30086,207.00629,106.21383,139.60162,185.87368,156.14022,209.69609
        ,90.36128,169.60823,97.66052,198.09048,65.52655,174.88725,214.64368
        ,137.93957
    };
    double y2[50] = {
        37.907311,180.827621,269.188404,188.385973,213.534848,213.860997
        ,140.896068,149.150317,184.048664,5.314890,240.332382,88.005440
        ,50.848480,211.425323,148.917861,212.530752,150.368739,112.340123
        ,7.367409,39.607834,45.064429,5.394230,176.623531,153.412619
        ,148.732398,88.662064,14.156160,3.841351,85.866524,98.224884
        ,174.691692,84.835332,81.432231,166.109418,120.217455,105.200469
        ,98.664827,38.200510,120.661322,299.502540,293.094608,17.254362
        ,190.781908,141.508701,67.785361,112.853815,6.877680,118.937008
        ,148.874783,242.670987
    };
    
    vector<double> vx1(x1, x1+20);
    vector<double> vy1(y1, y1+50);
    
    vector<double> vx2(x2, x2+50);
    vector<double> vy2(y2, y2+50);
    
    cout<<welchTest(vx1, vy1, criticalValuesAt80)<<endl;
    cout<<welchTest(vx2, vy2, criticalValuesAt99)<<endl;
    
    bool testRes1 = true;
    //significance level = 0.2
    testRes1 &= welchTest(vx1, vy1, criticalValuesAt80) == 1;
    //significance level = 0.1
    testRes1 &= welchTest(vx1, vy1, criticalValuesAt90) == 0;
    //significance level = 0.05
    testRes1 &= welchTest(vx1, vy1, criticalValuesAt95) == 0;
    //significance level = 0.025
    testRes1 &= welchTest(vx1, vy1, criticalValuesAt975) == 0;
    //significance level = 0.01
    testRes1 &= welchTest(vx1, vy1, criticalValuesAt99) == 0;

    bool testRes2 = true;
    //significance level = 0.2
    testRes2 &= welchTest(vx2, vy2, criticalValuesAt80) == 1;
    //significance level = 0.1
    testRes2 &= welchTest(vx2, vy2, criticalValuesAt90) == 1;
    //significance level = 0.05
    testRes2 &= welchTest(vx2, vy2, criticalValuesAt95) == 1;
    //significance level = 0.025
    testRes2 &= welchTest(vx2, vy2, criticalValuesAt975) == 1;
    //significance level = 0.01
    testRes2 &= welchTest(vx2, vy2, criticalValuesAt99) == 0;
    
    return testRes1 && testRes2;
}


}
