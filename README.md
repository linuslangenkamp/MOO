    ```textfile
                                          __  __   ___   ___
                                         |  \/  | / _ \ / _ \   MODELICA
                                         | |\/| || | | | | | |     / OPTIMIZER
                                         | |  | || |_| | |_| |    /
                                         |_|  |_| \___/ \___/   🐄 
                  
                  +-------------------------------------------------------------------------+
                  |          External Modeling Environment, e.g., OpenModelica (OM)         |
                  |               (Generates problem, requests optimization)                |
                  +-------------------------------------------------------------------------+
                       |
                       | (Problem definition, sparsity, callbacks for evaluations)
                       v
                  +-------------------------------------------------------------------------+
                  |   interfaces/openmodelica/                                              |
                  |   - gdop_problem.* (GDOP problem definition, etc.)                      |
                  |   - evaluations.* (Callbacks for objective, Jacobians, Hessians)        |
                  |   - main_opt.cpp (Entry point for OM-generated code)                    |
                  |   - sim_runtime_ext.* (External OM simulation runtime code)             |
                  |   - strategies.* (Special OM strategies / overrides)                    |
                  +------------------------+------------------------------------------------+
                                         |
                                         | (Optimization Request & Problem Details)
                                         v
                  +-------------------------------------------------------------------------+
                  |   nlp/instances/gdop/                                                   |
                  |   +-----------------------------------------------------------------+   |
                  |   |  gdop_orchestrator.* (The Brain for GDOP)                       |   |
                  |   |  - Decides *how* to solve GDOP (flow control, strategy choice)  |   |
                  |   |  - Manages overall optimization process (e.g., mesh refinement) |   |
                  |   +----+------------------+------------------+---------------------+    |
                  |        |                  |                  |                          |
                  |        v                  v                  v                          |
                  |   +-------------------+ +----------------------+ +--------------------+ |
                  |   | gdop_strategies.* | | gdop.*               | | problem.h          | |
                  |   | (Specific traits) | | (NLP implementation) | | - Base problem     | |
                  |   | - Mesh refinement | | - Defines GDOP - as  | | - Interface        | |
                  |   | - Initialization  | |   an NLP             | +--------------------+ |
                  |   | - Scaling logic   | | - Handles callbacks, |                        |
                  |   | - etc.            | |   evaluations        |                        |
                  |   +-------------------+-+----------------------+                        |
                  |                            |                                            |
                  |                            | (Uses generic NLP interface)               |
                  |                            v                                            |
                  |   +-----------------------------------------------------------------+   |
                  |   |  nlp/                                                           |   |
                  |   |  - nlp.h, nlp_solver.h (Generic NLP interface, structure)       |   |
                  |   |  - nlp_scaling.* (Generic scaling methods)                      |   |
                  |   |  +------------------------------------------------------------+ |   |
                  |   |  | solvers/                                                   | |   |
                  |   |  | - ipopt/ (ipopt_adapter.*, ipopt_solver.* Ipopt interface) | |   |
                  |   |  | - etc. WIP                                                 | |   |
                  |   |  | - nlp_solver_flags.* (Generic solver configuration flags)  | |   |
                  |   |  +------------------------------------------------------------+ |   |
                  |   +-----------------------------------------------------------------+   |
                  +-------------------------------------------------------------------------+
                                               |
                                               | (Base Operations & Helpers)
                                               v
                  +-------------------------------------------------------------------------+
                  |   base/                                                                 |
                  |   - linalg.*, mesh.*, trajectory.* (Core mathematical & numerical ops)  |
                  |   - collocation.*, util.h, log.h (Helper utilities)                     |
                  +-------------------------------------------------------------------------+
    ```
