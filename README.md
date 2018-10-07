# Optimizer
C++ solutions to optimization problems

## About
The goal is to build an Optimization library for C++ which can provide solutions for single variable and muli-variable optimization, constrained and unconstrained problems. The users should be able to select the algorithm to be used and also access any intermediate data which is obtained whilst running the algorithm. 

We want the contributors to learn how C++ libraries work and also learn some mathematics. Some minimal efforts are required for each good PR. Unlike other repositories hosted for Hacktoberfest, we want to create a fully utilizable product for developers and researchers over the months of October and beyond. If you wish to be a part of this we promise you that your efforts would bring over a good impact to everyone who will be utlizing this product. Best of Luck!

**This repository is also hosted as a part of the Hacktoberfest. Register yourself [here](https://hacktoberfest.digitalocean.com/) and submit 5 PRs to win a free T-shirt. This repository is for you to learn the github workflow, understand the community and get started with Open Source.**

## Understanding the Repository
- [src](Optimizer/src) contains all the code
- [Theory](Theory) contains all the algorithms we need to implement
- [Features](Features) has information about current features as well as features required. **New Contributors head here!**

## How to get this working?
**Currently there is no support for Windows.** For MAC/ Linux users :
- Fork and then clone this repository from your fork
- Install `Eigen` if you haven't already. Follow the instructions given [here](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)
  - Once you have downloaded the zip file from the link, extract it
  - Inside the extracted folder there will be a folder called `Eigen`
  - Cut and paste that into your `/usr/local/include`
  - If you have done `brew install Eigen` (or the equivalent on Linux), you can check your `usr/local/include` and open the     `eigen3` folder and move the `Eigen` folder here outside to `usr/local/include`(If this folder already has `Eigen` then     you may skip this step).
  - Go to the `Tests` folder and compile and run `test.cpp`. If there is no error then you have succesfully installed           `Eigen`, if not then you may get something like `Eigen/Dense: No such file or directory`. Revisit the previous steps and      see if you have screwed up somewhere.
- Inside [test.cpp](Tests/test.cpp) you can see the line `#include "../Optimizer/optimizer"` on top. This is a specific link to the Optmizer  header file. You can see that if you move [test.cpp](Tests/test.cpp) somewhere outside this folder, it will not work.
- To make sure that any `.cpp` file you use can include the Optmizer library by adding the line : `#include <Optimizer/optimizer>` you need to move the [Optmizer](Optimizer) folder in this repository into your `usr/local/include`. Once that is done check if it is working by moving [test.cpp](Tests/test.cpp) to some other arbitrary folder and changing this : `#include "../Optimizer/optimizer"` to this : `#include <Optimizer/optimizer>`
- In case there are any errors, then make sure you report them as issues. Also before you report them, do check if you are using `C++11`. Mac OS `g++` compilers use an older version of C++. You may have to compile the code like this : 
`g++ test.cpp -o test -std=c++11` 
 

## Key concepts
If you are pursuing engineering/ will be pursuing engineering/ have already pursued engineering, then this should be very easy for you to grasp. Without going too much into the mathematics you can easily code up the algorithms if you wish. But as a rule of thumb familiarity with these concepts would be helpful : 
- Basic C++
- Basic Calculus
- Linear Algebra
- Polynomials and Polynomial Systems
- Linear Programming
- Quadratic Programming

**DON'T BE INTIMIDATED! YOU DON'T HAVE TO STUDY ALL THIS!** Read on to get started ...

## How to Contribute
- Read this [README](README.md) completely.
- Checkout [Features](Features/README.md) and [Issues](https://github.com/sudz123/Optimizer/issues)
- Choose one which interests you, if it's already an Issue then add a comment claiming it. Else create a new issue.
- Read the respective section from [Theory](Theory) to understand how to get the code right
- Add the code and submit the PR by following the github workflow which you can find [here](https://egghead.io/lessons/javascript-how-to-fork-and-clone-a-github-repository)

## Your first PR!
Go and add your name to the [contributors.md](contributors.md) file.

## Important
- Kindly do look into the [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md)
- Follow the repository guidelines given in [RULES](Rules/README.md)
- Contact `sudarshan.kamath97@gmail.com`/ `diptanshu.agarwal@gmail.com ` for any issues(don't spam)
