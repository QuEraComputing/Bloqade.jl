# The Julia Programming Language

The Bloqade project is built in the [Julia](https://julialang.org) programming language. 
For those who are not familiar with Julia, here is a quick start for some basic Julia grammar,
and a guide for learning more about Julia and advanced usage.

## Why Julia?

- **Fast**. As you might have heard, Julia is very fast; there are various benchmarks online.
    It can be used to even write Basic Linear Algebra Subroutine (BLAS) to reach performance that is on par with the
    manually optimized assembly with C (check [Octavian](https://github.com/JuliaLinearAlgebra/Octavian.jl)).
- **Generic**. The language itself and its ecosystem are built to be generic, and the compiler can specialize
    on generic methods automatically. Thus, you will find that a lot things can be combined easily, and they will
    **just work**, e.g. plugging in the `Measurement` number from 
    [Measurement.jl](https://github.com/JuliaPhysics/Measurements.jl) into your ODE solver, you will get error propagation **automatically**; plugging in Tropical numbers or in general a semi-ring algebra into a tensor-network contraction function, you can
    [solve optimization problems with tensor networks](https://github.com/QuEraComputing/GenericTensorNetworks.jl), and so on.
- **Differentiable**. The language is differentiable, which means you can calculate the derivatives
    using an automatic differentiation (AD) engine on the whole language. The AD ecosystem in Julia is very well developed and supported. 
    The current stable AD engine
    is powered by [Zygote](https://arxiv.org/abs/1907.07587) and the next generation AD engine includes
    [Diffractor (check the video talk on ACM SIGPLAN)](https://www.youtube.com/watch?v=mQnSRfseu0c) and
    [Enzyme](https://enzyme.mit.edu).
- **Extensible**. The language is designed to be compiler friendly. It supports staged programming
    as well as compiler plugins. This makes supporting new hardware much easier. As a result, Julia
    can conveniently support multiple different hardware, such as [CUDA](https://github.com/JuliaGPU/CUDA.jl),
    [oneAPI](https://github.com/JuliaGPU/oneAPI.jl), [TPU](https://github.com/JuliaTPU/XLA.jl), and potentially quantum computers in the future.
- **Easy**. With all these powerful features, the language itself is still rather easy to learn. Let's go to
    the [quick start](#Quick-Start) section to skim through the basic syntax.


!!! info

    Multi-stage programming (MSP) is a variety of metaprogramming in which compilation is divided into a series of intermediate phases, allowing type-safe run-time code generation. Statically defined types are used to verify that dynamically constructed types are valid and do not violate the type system. 

    -- Wikipedia

## Quick Start

### Variables and Some Basic Types
In Julia, you can define a variable similar to how you define it in Python. 
For example, you can define a `x` using `=` (assignment):

```@repl quick-start
x = 1
```

Every variable has a type. You can check it using `typeof`:

```@repl quick-start
typeof(x)
```

By default, Julia displays the output of the last operation. You can suppress the output by adding `;` (a semicolon) at the end.

### Functions
In Julia, you can also define short-form, one-line functions using `=` (assignment) similar to how you write things mathematically.

```@repl quick-start
f(x) = 2x
```

Typing the function's name gives information about the function. To call it, we must use parentheses:

```@repl quick-start
f
f(2)
```

For longer functions, we use the following syntax with the `function` keyword and `end`:

```@repl quick-start
function g(x, y)
	z = x + y
	return z^2
end
```

### Control Flows
In Julia, there are `for`, `if` and `while` control flows. For example, the `for` loop looks like:

```@repl quick-start
s = 0
for i in 1:10
    s += 1
end
```

we can now check the value of `s` by typing it again:

```@repl quick-start
s
```

Here, `1:10` is a **range** representing the numbers from 1 to 10:

```@repl quick-start
typeof(1:10)
```

the `if else` statement looks like the following:

```@repl quick-start
if s < 10
	# do something
elseif 10 < s < 13
	# do something
else
	# do something
end
```

### Matrix and Array
Julia carries its own `Array` type. If you use Python, it is similar to `numpy.array` in Python except:
1. index starts from 1,
2. the multi-dimensional index is column-wise.
You can also use list comprehension:

```@repl quick-start
[i for i in 1:10]
```

It works for multi-dimensional case too:

```@repl quick-start
[(i, j) for i in 1:10, j in 1:5]
```

Most functions involving matrices and arrays follow the same convention as `numpy` or `MATLAB`. For example, you can create a random matrix using:

```@repl quick-start
rand(5, 5)
```

If you have questions about using a function, you can always type the question mark `?` in your REPL following the function name:
```julia
julia> ?rand
```

### Package Manager & Environments
Julia carries its own package manager. You can use it as a normal package:
```julia
julia> using Pkg
```
To install a package, you can use:
```julia
julia> Pkg.add("Bloqade")
```
To remove a package, you can use:
```julia
julia> Pkg.rm("Bloqade")
```
All Julia programs **run inside an environment**. The default is the global environment. It is usually recommended to run your notebook in a local environment, so you won't hit **any version conflicts** between different packages. 

## Resources

For more resources, check the official website [julialang.org/learning](https://julialang.org/learning):

```@raw html
<style>
  iframe {
    width: 1px;
    min-width: 100%;
    min-height: 1000px;
  }
</style>
<iframe id="myIframe" src="https://julialang.org/learning/"></iframe>
<script>
  iFrameResize({ log: true }, '#myIframe')
</script>
```