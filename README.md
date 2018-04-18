```
//     _____________________________________      _____   |    
//     ___/ __ \__/ __ \_/ __ \__/ __ \__/ /________/ /   |
//     __/ /_/ /_/ /_/ // / / /_/ /_/ /_/ __ \/ _ \/ __/  |
//     _/ ____/_/ _, _// /_/ /_/ ____/_/ / / /  __/ /_    |
//     /_/     /_/ |_| \____/ /_/     /_/ /_/\___/\__/    |
//---------------------------------------------------------
```
[![DOI](https://zenodo.org/badge/68395684.svg)](https://zenodo.org/badge/latestdoi/68395684)

 PROPhet was ceated by      
 Brian Kolb and Levi Lentz     
  in the [group of Alexie Kolpak](http://kolpak.mit.edu/) at [![MIT](https://github.com/biklooost/PROPhet/blob/gh-pages/images/MIT-logo-red-gray-54x28.png)](http://web.mit.edu/)                                   

  See the [documentation](https://github.com/biklooost/PROPhet/blob/master/doc/PROPhet_documentation.pdf) for
  a full description of PROPhet.

  Please see the [LICENSE](https://github.com/biklooost/PROPhet/blob/master/LICENSE) file for license information
  
  If used for published work, please cite the following:
  
  [Scientific Reports 7, Article number: 1192 (2017)
   doi:10.1038/s41598-017-01251-z](https://www.nature.com/articles/s41598-017-01251-z)
                                                           
---


### [Description of PROPhet](https://biklooost.github.io/PROPhet/)
---
PROPhet (short for PROPerty Prophet) couples neural networks with first-principles physics and chemistry codes to 
allow sophisticated prediction of material properties. In general, PROPhet is used to
find mappings between a set of material or system properties and other properties.  Some
specific uses of PROPhet are to: 
* Find a mapping between atomic configuration and other properties, including the total energy, creating an analytical potential, which can be used for molecular dynamics with the [LAMMPS](http://lammps.sandia.gov/) code
* Construct density functionals for exchange-correlation energy, kinetic energy, or just about
anything else
* Find a mapping between a set of descriptors (an arbitrary combination of material properties) and other properties

See the [PROPhet project page](https://biklooost.github.io/PROPhet/) for more details.


### Interfaces with other codes
---
At the moment, PROPhet can couple automatically to the first-principles codes
- [VASP](https://www.vasp.at/)
- [Quantum Espresso](http://www.quantum-espresso.org/)
- [FHI-Aims](https://aimsclub.fhi-berlin.mpg.de/)

meaning it can extract many common properties directly from the output files of these codes, without user interaction.  Interfaces to other codes can be easily added.

In addition, potentials created in PROPhet can be used for molecular dynamics runs with the freely-available [LAMMPS](http://lammps.sandia.gov/) MD code.  See the [documentation](https://github.com/biklooost/PROPhet/blob/master/doc/PROPhet_documentation.pdf#page=16) for more details.



### Compiling Information
---
Compilation follows the standard linux paradigm:
```
./configure [options]
make
make install
```

If you want to use PROPhet potentials in the [LAMMPS](http://lammps.sandia.gov/) MD code, 
you should execute:

```
./configure [options] --enable-lammps=LAMMPS_DIR
make 
make install
```

where LAMMPS_DIR is the directory with the LAMMPS source.
This will make the lammps library of PROPhet, and attempt 
to insert it into the LAMMPS package system.  The LAMMPS
code must be relinked after completion to link in the 
PROPhet library.  If automatic instalation fails, the 
library can be inserted into LAMMPS by following the 
instructions given in the [documentation](https://github.com/biklooost/PROPhet/blob/master/doc/PROPhet_documentation.pdf#page=6).

### Usage Instructions
---
Usage instructions including a tutorial can be found in the [documentation](https://github.com/biklooost/PROPhet/blob/master/doc/PROPhet_documentation.pdf#page=8),
and by looking in the [doc/tutorial](https://github.com/biklooost/PROPhet/tree/master/doc/tutorial) directory. A fully
annotated input_file is also included in the 
[doc](https://github.com/biklooost/PROPhet/tree/master/doc) directory.

