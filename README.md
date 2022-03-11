# StochasticGrowthRate

RUNNING SENSTOCH.R NOTES

TO DO:
--check P1 matrix in Pascarella and Horvitz vs cmat from code

For Mac OSX:
Download Xcode
Download command line tools (if not already included in Xcode installation)

Download R packages:
lpSolve
expm
inline
Rcpp
RcppArmadillo*

*The more recent versions of the RcppArmadillo R package do not have the correct eig_gen function (necessary for the C++ senstoch.cpp code). To get around this, install RcppArmadillo version 4 from the tar.gz.


Keep senstoch.cpp and senstoch.R in the same folder--senstoch.R sources senstoch.cpp

After sourcing senstoch.R, DR's machine has error messages, but still runs fine. Test using function e.g. statdist(P1). Should output a vector (stationary distribution of P1 matrix)

Note from David 9/14/13:
Here is the first complete version of the stochastic growth software. The past few weeks I ran some final tests, where I actually checked the output of the sensitivity routines against the differences you get from small changes to the matrices. They didn't match up at all, and for a while I had convinced myself that our published results were wrong. But, hooray for rigorous mathematics, when I took the results apart again, I realised that the methods really were solid. But programming is even less forgiving of sloppiness than mathematics, so it took a few weeks to track down a mismatched index here, a flipped sign there. But now it seems to be producing consistent results on several different examples. And it's pretty fast.

I haven't gotten to write up anything like formal documentation, but here are some basic descriptions of how to use the program.

First of all, installation isn't as trivial as one might hope. The program is, at this point, basically in C++, with just a thin veneer of R to keep you from going mad. Until this gets turned into a package (that can be downloaded as a binary already compiled) that means that your computer needs to be able to compile C code.

If you're running Linux, you may not need to do anything special. But I haven't checked.

If you're running Mac OSX (which is what I have), you need to first download Xcode (free from the app store, or directly from the Apple developer web site, or on the installation DVD that comes with every new computer). Then you probably need the command-line tools (which used to be part of Xcode, but now need to be downloaded separately. Once Xcode is installed, run it, and go to the preferences pane, and then to Downloads, and click on Install by command-line tools. Then you need to download and install gfortran from http://cran.r-project.org/bin/macosx/tools/

You need to download and install the following standard R packages from CRAN: lpSolve, expm, inline, Rcpp, RcppArmadillo.

When you're finished with that, put the two files into the same directory, and source senstoch.R. It will source senstoch.cpp and compile it, which will take a bit of time (something like 15 seconds on my machine). For some reason, on one of my machines it gives a series of warnings of various directories not found, but this doesn't seem to affect the functioning. You can try typing in a command like
statdist(P1),
which outputs the stationary distribution of the stochastic matrix P1. If that outputs a vector and not an error message, the compilation has probably gone alright.

After that, you're ready to go. The only commands you really need are

SGmean(X, P, NU):  Computes the stochastic growth rate
sensitivity1(X, P, NU): Computes all sensitivities with respect to changes in the projection matrices
sensitivity2(X, P, N, NU): Computes all sensitivities with respect to changes in the environment transitions

The inputs are X= list of projection matrices. All need to be the same (square) size.
P= transition matrix on environments. This needs to be K x K, where K is the length of the list X.
NU= # of sample vectors taken to represent the stationary distribution of the age distribution. (Left and right vectors are, by default, sampled in the same numbers.)
N = # of samples used in estimating zeta. This is a total number of samples, spread over all possible coupling times and states, and no matter what you choose it's going to take 25 each at a minimum. It chooses the maximum coupling time to get the estimated bias down below a given threshold. When working with these 8x8 matrices it goes for a maximum coupling time of about 100 (this is something I need to think about optimising some more), so it should be a minimum of about 20000. (In the example below I only entered 10000, so it was probably irrelevant. Running it with N=100000 led to serious machine slowdown.)

Sensitivity calculation time (and memory requirements) increase with the square of NU, so it shouldn't be much bigger than 1000. One problem with running C code is that it seems to have no natural restrictions on how much of your computer resources it can take over. (I need to look into this.) So if you make your N's too big, you can end up with R running forever, and nothing else functioning on your computer either (at least, on the Mac). A couple of times I've needed to restart my computer because of this problem, and numerous times I needed to shut down R.

You should probably start with conservative choices of N and NU (NU=100 for sensitivities, up to 1000000 for SGmean, and N=10000). If the variances are too large, you can either choose a larger N and NU, or you can just run the same thing multiple times and average the results. This will reduce the variance, and leave the bias unchanged.

There are some other inputs to these functions that you can probably just leave at the default settings to begin with.

All numbers are output with a mean, a variance, and a bias. The mean is the point estimate for the quantity. A 99.7% confidence interval for the quantity in question will then be

Mean +/- (Bias + 3*sqrt(Variance))

If you want a different sort of confidence interval you replace 3 by something else.

SGmean outputs just three numbers.

sensitivity1 outputs a list, whose length is the number of matrices in X. For each d x d matrix in X there is a Mean matrix, a Variance matrix, and a Bias matrix, each one also a d x d matrix, showing the sensitivity of stochastic growth to changes in the corresponding entry of the X matrix.

sensitivity2 outputs a list of three K x K matrices, corresponding to Mean, Variance, and Bias. The Mean matrix shows the sensitivity of stochastic growth to changes in the corresponding entry of P. Because only differences within a row have any meaning I have normalised to make one entry in every row 0, and all other entries positive.

Let me know how it goes. Probably there will be some immediate bugs that will come up as soon as you try to run it on a different machine.
[see attached file: senstoch.R]
[see attached file: senstoch.cpp]

I programmed in the matrices from the paper with the seven 8x8 matrices at the end of the file senstoch.R. The following computation took about 10 seconds on my several-year-old laptop.

> s2=sensitivity2(X,P1,10000,500)
> s2
$mean
     [,1]       [,2]       [,3]      [,4]      [,5]      [,6]      [,7]
[1,]    0 0.39819106 1.06440920 1.5455268 2.2138087 3.2228336 4.4775992
[2,]    0 0.04138872 0.09680221 0.1559391 0.2445112 0.3597096 0.5081297
[3,]    0 0.03713765 0.08087182 0.1190507 0.2104891 0.3149138 0.4445068
[4,]    0 0.04017335 0.08115001 0.1241095 0.1807485 0.2791917 0.3823390
[5,]    0 0.02553069 0.06972826 0.1299994 0.2054757 0.3118477 0.4336626
[6,]    0 0.05504393 0.12053052 0.2144009 0.3301829 0.4886798 0.6775851
[7,]    0 0.03852828 0.08686189 0.1794192 0.2953095 0.4424452 0.6160436

$variance
     [,1]         [,2]         [,3]         [,4]         [,5]         [,6]         [,7]
[1,]    0 6.294626e-04 5.779874e-04 6.108461e-04 6.370242e-04 6.124284e-04 6.677312e-04
[2,]    0 1.655878e-05 1.984939e-05 1.887580e-05 1.853706e-05 1.730501e-05 1.944866e-05
[3,]    0 1.331835e-05 1.506914e-05 1.483163e-05 1.388791e-05 1.382952e-05 1.398683e-05
[4,]    0 1.069623e-05 1.137104e-05 1.033353e-05 1.037369e-05 1.097223e-05 1.174663e-05
[5,]    0 1.328988e-05 1.394504e-05 1.410879e-05 1.439337e-05 1.407199e-05 1.556999e-05
[6,]    0 1.846955e-05 2.117830e-05 1.851712e-05 1.877111e-05 1.911898e-05 1.775998e-05
[7,]    0 2.943529e-05 3.140280e-05 3.047464e-05 2.813832e-05 2.812648e-05 2.712696e-05

$bias
     [,1]         [,2]         [,3]         [,4]         [,5]         [,6]         [,7]
[1,]    0 3.058810e-04 2.965343e-04 3.064996e-04 3.266349e-04 3.350635e-04 3.771617e-04
[2,]    0 3.421539e-05 3.356554e-05 3.488219e-05 3.806640e-05 4.089013e-05 4.289307e-05
[3,]    0 2.981940e-05 3.069929e-05 3.142247e-05 3.347204e-05 3.495817e-05 3.807631e-05
[4,]    0 2.431428e-05 2.397463e-05 2.623044e-05 2.710588e-05 2.883903e-05 3.202155e-05
[5,]    0 2.950036e-05 3.019744e-05 3.056974e-05 3.102973e-05 3.373306e-05 3.917898e-05
[6,]    0 4.583045e-05 4.595810e-05 4.778320e-05 4.914067e-05 5.378836e-05 5.829750e-05
[7,]    0 4.444509e-05 4.453496e-05 4.706111e-05 4.807976e-05 5.138786e-05 5.719820e-05
