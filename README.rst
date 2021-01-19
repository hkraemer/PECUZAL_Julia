.. image:: https://img.shields.io/badge/docs-dev-blue.svg
    :target: https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/unified/

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4449785.svg
   :target: https://doi.org/10.5281/zenodo.4449785

PECUZAL Julia
=============
This code base is using the Julia Language and `DrWatson <https://juliadynamics.github.io/DrWatson.jl/stable/>`_
to make a reproducible scientific project, authored by K.Hauke Kraemer and
George Datseris. It contains all the source code for producing the article [kraemer2020]_.


To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git repo and may need to be downloaded independently.
1. Open a Julia console and do:

   .. code-block:: julia

       julia> cd("path/to/this/project")
       julia> using Pkg; Pkg.activate(".")
       julia> Pkg.instantiate()

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

Documentation
=============
This repository is not intended to make a self-contained Julia package, but rather
give readers of [kraemer2020]_ the opportunity to get full access to any source
code.
The PECUZAL method is incorporated and maintained in the
`DynamicalSystems.jl <https://juliadynamics.github.io/DynamicalSystems.jl/dev/>`_-Ecosystem,
specifically in the `DelayEmbeddings.jl <https://github.com/JuliaDynamics/DelayEmbeddings.jl>`_
package. `Here <https://juliadynamics.github.io/DynamicalSystems.jl/latest/embedding/unified/>`_
the reader can find a full documentation and some basic example illustrating the usage of the PECUZAL method.

Citing and reference
====================
If you enjoy this tool and find it valuable for your research please cite

.. [kraemer2020] Kraemer et al., "A unified and automated approach to attractor reconstruction",  `arXiv:2011.07040 [physics.data-an] <https://arxiv.org/abs/2011.07040>`_, 2020.

or as BiBTeX-entry:

::

    @misc{kraemer2020,
    title={A unified and automated approach to attractor reconstruction},
    author={K. H. Kraemer and G. Datseris and J. Kurths and I. Z. Kiss and J. L. Ocampo-Espindola and N. Marwan},
    year={2020},
    eprint={2011.07040},
    archivePrefix={arXiv},
    primaryClass={physics.data-an}
    url={https://arxiv.org/abs/2011.07040}
    }


Licence
=======
This is program is free software and runs under `MIT Licence <https://opensource.org/licenses/MIT>`_.
