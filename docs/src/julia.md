# The Julia Programming Language

The Bloqade project is built in pure [Julia](https://julialang.org) programming language. Thus
if you are not familiar with Julia, here is a quick start for basic Julia gramar,
and a guide for learning more detailed and advanced Julia.

## Why Julia?

- **Fast**, as you might have heard about it, Julia is very fast, there are various benchmarks online.
    It can even be used to write Basic Linear Algebra Subroutine (BLAS) to reach performance on par with
    manually optimized assembly with C (check [Octavian](https://github.com/JuliaLinearAlgebra/Octavian.jl)).
- **Generic**, the language itself and its ecosystem are built to be generic, and the compiler can specialize
    on generic methods automatically, thus you will find a lot things can be combined easily, and they will
    **just work**, e.g plugin the `Measurement` number from 
    [Measurement.jl](https://github.com/JuliaPhysics/Measurements.jl) into your ODE solver, you will get error propagation **just work**, plugin Tropical number into tensor contraction function, you can
    [solve optimization problems with tensor networks](https://arxiv.org/abs/2008.06888), and so on.
- **Differentiable**, the language is differentiable, that means you can calculate the derivatives
    using an automatic differentiation engine on the whole language. The current stable AD engine
    is powered by [Zygote](https://arxiv.org/abs/1907.07587), the next generation AD engine includes
    [Diffractor (check the video talk on ACM SIGPLAN)](https://www.youtube.com/watch?v=mQnSRfseu0c),
    [Enzyme](https://enzyme.mit.edu).
- **Extensible**, the language is designed to be compiler friendly, it supports staged programming
    as well as compiler plugins. This makes supporting new hardware much easier. As a result, Julia
    can support multiple different hardware, such as [CUDA](https://github.com/JuliaGPU/CUDA.jl),
    [oneAPI](https://github.com/JuliaGPU/oneAPI.jl), [TPU](https://github.com/JuliaTPU/XLA.jl) and so on.
- **Easy**, with all these features, yet the language itself stays rather easy to learn. Let's go to
    the [quick start](#Quick-Start) section to skim the syntax.


!!! info

    Multi-stage programming (MSP) is a variety of metaprogramming in which compilation is divided into a series of intermediate phases, allowing typesafe run-time code generation. Statically defined types are used to verify that dynamically constructed types are valid and do not violate the type system. -- Wikipedia

## Quick Start

### Variables and Some Basic Types
In Julia, you can define a variable similar to how you define it in Python, e.g we can define a `x` using `=` (assignment)

```@repl quick-start
x = 1
```

every variable has a type, you can check it using `typeof`

```@repl quick-start
typeof(x)
```

By default Julia displays the output of the last operation. (You can suppress the output by adding `;` (a semicolon) at the end.)

### Functions
In Julia, you can also define short-form, one-line functions using `=` (assignment) similar to how you write things mathematically.

```@repl quick-start
f(x) = 2x
```

Typing the function's name gives information about the function. To call it we must use parentheses:

```@repl quick-start
f
f(2)
```

For longer functions we use the following syntax with the `function` keyword and `end`:

```@repl quick-start
function g(x, y)
	z = x + y
	return z^2
end
```

### Control Flows
In Julia, there are `for`, `if` and `while`, they look like the following

```@repl quick-start
s = 0
for i in 1:10
    s += 1
end
```

we can now check the value of `s` by typing it again

```@repl quick-start
s
```

Here, `1:10` is a **range** representing the numbers from 1 to 10:

```@repl quick-start
typeof(1:10)
```

the if else statement looks like the following

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
Julia carries its own `Array` type, if you use Python, it is similar to `numpy.array` in Python except:
1. index starts from 1
2. the multi-dimensional index is column-wise
You can also have list comprehension:

```@repl quick-start
[i for i in 1:10]
```

it works for multi-dimensional case too:

```@repl quick-start
[(i, j) for i in 1:10, j in 1:5]
```

most functions follow the same convention as numpy or MATLAB, e.g you can create a random matrix using:

```@repl quick-start
rand(5, 5)
```

if you have question about using a function, you can always type question mark `?` in your REPL following the function name
```julia
julia> ?rand
```

### Package Manager & Environments
Julia carries its own package manager, you can use it as a normal package:
```julia
julia> using Pkg
```
to install a pacakge, you can use
```julia
julia> Pkg.add("Yao")
```
to remove a pacakge, you can use
```julia
julia> Pkg.rm("Yao")
```
All Julia program **runs inside an environment**, it is the global environment by default. It is usually recommended to run your notebook in a local environment, so we won't hit **any version conflicts** between different packages. 

## Resources

For more resources just check the official website [julialang.org/learning](https://julialang.org/learning)

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