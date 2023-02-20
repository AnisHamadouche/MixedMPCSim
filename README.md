## MixedMPCSim Benchmark

MixedMPCSim is an open-source MATLAB benchmark for reduced precision MPC simulation. MixedMPCSim is based on reguralized optimization formulation of MPC. You can choose between two numeric solvers: the alternating direction method of multipliers (ADMM) or the Proximal-gradient Descent (PGD) under different floating-point and fixed-point machine representations.  

**Current version:** 1.0.0

**Release notes:** 

* The current version only implements PGD and ADMM for the <img src="https://latex.codecogs.com/svg.latex?\Large&space;\ell_1" title="\ell_1" /> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;\ell_0" title="\ell_0" /> regularized LASSO problem.

* Other first-order solvers and other problems will be added in future releases. 

## Contents
* [Description](#Description)
* [Quick start](#QuickStart)
* [How to cite](#References)
* [Contact us](#Contacts)
* [Licence](#Licence)


## Description<a name="Description"></a>

MixedMPCSim Benchmark implements PG and ADMM to solve the MPC problem

<img src="https://latex.codecogs.com/svg.latex?\Large&space;(1)\quad\quad\text{minimize}\quad\quad\frac{1}{2}\|Au-b\|_2^2+\|u\|_1" title="\Large \text{minimize}\quad\quad\frac{1}{2}\|Au-b\|_2^2 + \|u\|_p" />

where p = 1 or p = 0. MixedMPCSim Benchmark offers a choice to solve problem (1) using different custom data types.

## Requirement<a name="Requirement"></a>

* CVX
* Fixed-point Designer

## Quick start<a name="QuickStart"></a>

* Run
 	`` >> main ``

with custom problem data (add model, constraints and signals to ./models) or use the default spacecraft problem data to compare the output and control/differential control signal values under different machine representations. 

* The default benchmark runs PG and ADMM under 'double precision', 'single precision', '12 bits fixed-point' and  '16 bits fixed-point' representations. Enter custom data casting precision in "frmt.data" and solver precision in "frmt.dt". 

* Choose the problem via "para.problem" which can be set to either "l0" or "l1".
* Choose the solver via "para.solver"  which can be set to either "pgd" or "admm". 

To add custom data types add a case statement with custom data type name. To invoke a specific type within another function use  
	`` >> T = mixedTypes('data type'); ``
then use casting as follows:
	`` >> x = cast(x0, 'like', T.x) ``

	
**NOTE:** _this is a research code, and is under active development. You may find 
some undocumented inputs and options that are being used for development 
purposes, in the hope that they will become part of the "official" release. If 
you have any suggestions for improvement, or find any bugs, feel free to [contact us](#Contacts)!_


## How to cite<a name="References"></a>

If you find MixedMPCSIM Benchmark useful, please cite the following paper as appropriate:

```
@inproceedings{wu2020approximate,
  title={Approximate lasso model predictive control for resource constrained systems},
  author={Wu, Yun and Mota, Jo{\~a}o FC and Wallace, Andrew M},
  booktitle={2020 Sensor Signal Processing for Defence Conference (SSPD)},
  pages={1--5},
  year={2020},
  organization={IEEE}
}
	
```

## Contact us<a name="Contacts"></a>
To contact us about MixedMPCSim Benchmark, suggest improvements and report bugs, email either [Anis Hamadouche] (mailto:ah225@hw.ac.uk?Subject=MixedMPCSim).


## Licence<a name="Licence"></a>

MixedMPCSIM Benchmark is free software; you can redistribute it and/or modify it under the terms 
of the [GNU Lesser General Public Licence (LGPL)](https://www.gnu.org/licenses/lgpl-3.0.en.html) as published by the Free Software
Foundation; either version 3 of the Licence, or (at your option) any later version.

MixedMPCSIM Benchmark  is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. [See the GNU General Public License for more details](https://www.gnu.org/licenses/gpl-3.0.en.html).

You should have received a copy of the GNU Lesser General Public License along 
with MixedMPCSIM Benchmark; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
