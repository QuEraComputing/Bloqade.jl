# [Installation](@id install)

You can copy the following line to your Julia REPL
to install the latest stable version of this package:

```julia
pkg> add Bloqade
```

## Low-latency Usage of Bloqade Component Packages

The Bloqade project contains multiple packages. For development on top of the functionality,
(especially for those who do not need the ODE solvers), we recommend you to use the corresponding
component packages. The following is a list of component packages and what they do (WIP = work-in-progress)

- BloqadeExpr: Expressions and API definitions for Bloqade.
- BloqadeKrylov: Krylov-subspace based emulation.
- BloqadeLattices: objects, functions for lattices.
- BloqadeMIS: tools for working with maximum-independent sets in Rydberg system.
- BloqadeODE: ODE-based emulation.
- BloqadePython: WIP, python wrapper for the Bloqade package.
- BloqadeQMC: WIP, Stochastic Series Expansion for Rydberg system.
- BloqadeSchema: WIP, the schema for creating a task for Bloqade and QuEra machine.
- BloqadeWaveforms: the waveform objects.
- YaoSubspaceArrayReg: register object and functions in a subspace.

All the non-WIP packages are registered in the General registry. Thus, you can add them
as your dependency by directly running `pkg> add <component package>` in your Julia REPL.

## Try the Latest Version of Bloqade

Some users may want to try the latest version of Bloqade
for bug fixes, new features, etc. One can use `git` to clone the
repo to try the latest version of the entire package. This
requires one to setup the local project environment via `dev`.
Please refer to the page [Contributing to Bloqade](@ref) for more information.

If you only want to try the latest version of a specific
Bloqade package, just add `#master` after the package name, e.g.:

```julia
pkg> add BloqadeExpr#master
```

## Using Bloqade with AWS EC2 services

Bloqade is an open source solution for Hamiltonian simulation and can be deployed on any personal computer. Yet, some users might benefit from the extra performance offered by large computational resources from different providers. To address that, Bloqade is also available at the Amazon Web Services (AWS) Marketplace, and prepared to run on AWS EC2 instances via Amazon Machine Images. To deploy Bloqade on AWS EC2 instances, follow the 7-step process below.

- disclaimer 1: deploying Bloqade on AWS EC2 instances will incur a cost on the user that will depend on the AWS resources utilized.  
- disclaimer 2: support on deploying Bloqade on AWS can be obtained via AWS Support. This is a one-on-one support channel that is staffed 24x7x365 with experienced support engineers. To learn more, follow [this link](https://aws.amazon.com/premiumsupport/).  

### Step-0: Create an EC2 instance on AWS and check location

For the general guidelines on launching EC2 instances, check the [AWS EC2 tutorial](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html).
Bloqade can be run, technically, from any location. Still, its image is installed in servers held in N. Virginia. To get started as easily as possible, ensure to set the  AWS region location to N. Virginia(US-EAST-1)

<img width="745" alt="step0" src="https://user-images.githubusercontent.com/99290010/188176960-0fc2e132-d15c-4c64-90cd-99790a59b6ec.png">

### Step-1: Create an EC2 instance on AWS and check location

Now to really get started. On your AWS account portal, type EC2 on the search bar and access the EC2 service

<img width="654" alt="step1" src="https://user-images.githubusercontent.com/99290010/188177392-3681ff55-f418-44bd-9569-ee78ad9a272d.png">

### Step-2: Start the launch

Find the "launch instance" button to start running. You may check all the running instances clicking on the pointed button below.

<img width="411" alt="step2" src="https://user-images.githubusercontent.com/99290010/188177745-5c17d7f5-909d-4c72-846d-c724523df313.png">

### Step-3: Name your instance

Type a convenient name for your use...

<img width="771" alt="step3" src="https://user-images.githubusercontent.com/99290010/188178335-69a5bef0-2291-4989-9954-2686d52270ba.png">

### Step-4: Choose an image

...and choose Bloqade as an image. This will pre-install Julia, Bloqade, and Yao (for general quantum gates) on the cluster

<img width="495" alt="step4" src="https://user-images.githubusercontent.com/99290010/188178496-60d68d90-a5f1-42cc-81f5-2cc903496266.png">

### Step-5: Tune your settings

Select the EC2 instance type. Note that your charge rate will be function of these; large RAM charges more. For simple uses, we recommend m2.xlarge as a basic choice.

<img width="568" alt="step5" src="https://user-images.githubusercontent.com/99290010/188178728-f92e959e-a910-4e84-9bb6-b3872ade6cfc.png">

### Step-6: Tune more settings

Select your security group. This depends on your personal, company security operations, or AWS best practices.

<img width="735" alt="step6" src="https://user-images.githubusercontent.com/99290010/188178844-5cc2f714-fe14-4fac-aef9-9432c24ae0dc.png">

### Step-7: Blast off!

Launch your instance and Bloqade away!

<img width="652" alt="step7" src="https://user-images.githubusercontent.com/99290010/188178941-484e5c3d-a457-49c1-8944-e4417e8f820f.png">

### Final comments:

#### Bloqade Optimized AMIs from AWS Marketplace
The Bloqade team has implemented 2 dedicated AMIs that can be acquired using AWS Marketplace

##### Bloqade AMI with Julia (base image)

The Bloqade AMI contains:
- The latest Julia installation using `juliaup add release` which setup the latest Julia release. For more information please refer to [Juliaup](https://github.com/JuliaLang/juliaup) 
- The latest version of Bloqade 

##### Bloqade Optimized Deep Learning AMI with CUDA

In addition to the content of our base image, this image further contains
- NVIDIA CUDA, cuDNN, NCCL, GPU Drivers, Intel MKL-DNN, Docker, NVIDIA-Docker and EFA support
- Block devices


#### SSH Access

- On step 7, a key-pair login may be created for ssh access, if the user does not posses one previously.
- In this case, it is recommended to access Key pair (login) and click “Create key-pair”. Following the instructions, one may download a ````.pem```` file which can be added to the user’s ````~/.ssh```` folder. 
- Run the following command chmod ````400 <yourkey>````
- Following this process, a Defaul SSH protocol for sign in and security can be used.)
    * Do so by accessing ````vim config```` on your shell and typing the following
 <img width="612" alt="image" src="https://user-images.githubusercontent.com/99290010/188179655-dc6e8cbe-bdd0-462b-8d6a-dbd9cfccda7f.png">
You can find your DNS name  by clicking your instance on the list of running instances in the dashboard, clicking “Connect” and copying the “Public DNS” in the SSH Client tab

Then exit vim and type ````SSH AWS````

#### Other packages
- Bloqade runs in Julia but if you need to install extra Python packages for your uses, we recommend installing ````conda```` to manage virtual environments and install packages via ````conda install <package>````
- To install it via command line, follow 

````
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
````
You will need to exit out (````Ctrl-D````) of the AWS SSH ran re-log back in (```` ssh AWS````) to get ````conda```` activated.

<!---```
#/dev/sda1=snap-03d72fbeb983a4663:60:true:gp2
#/dev/sdb=ephemeral0
#/dev/sdc=ephemeral1
``` --->

  
## Build System Image to Accelerate Start-up Time

Since Bloqade is a large package, its loading time
and time-to-first-simulation can be very long.
You can build system images to save all the compilation
results in a binary to accelerate its loading/compilation
time. This is useful when you have lots of interactive
programming needs with Bloqade.

To build a system image for your environment, please use
the [PackageCompiler](https://julialang.github.io/PackageCompiler.jl/dev/)
or use the Julia VSCode plugin's [build system image feature](https://www.julia-vscode.org/docs/stable/userguide/compilesysimage/)


## Developing Bloqade

When developing Bloqade, one will need to setup a local environment
that contains all the local changes. To work with the Bloqade repo,
first you need to clone this repo

```sh
# clone this repo
git clone https://github.com/QuEraComputing/Bloqade.jl.git Bloqade
# go into the directory
cd Bloqade
# dev the corresponding environment
.ci/run dev
```

How does this work? The `.ci/run dev` command actually calls the `Pkg.develop`
command from Julia's package manager. Because we want to use the local
changes of the package, one will need to `dev` the corresponding package to 
make the changes happen in your current environment, e.g one will need to `dev` 
the `lib/BloqadeExpr` package to apply changes in `BloqadeExpr` module.

We also provide a convenient tool to setup this more automatically by
looking up dependencies in `lib` in one's `Project.toml` file,

```sh
.ci/run dev <path/to/your/environment>
```

will `dev` all the Bloqade dependencies in your environment.

See also [Modifying A Dependency](https://pkgdocs.julialang.org/v1/getting-started/#Modifying-A-Dependency)
for more detailed explainations.
