# Rydberg Emulator

The QuEra Rydberg System Emulator.

## Installation

If you have already added the private QuEra registry [Miskatonic](https://github.com/happy-Diode/Miskatonic), just type the following in your REPL pkg mode

```julia
(pkg)> add EaRydCore
```

To install the emulator, you need to add QuEra registry [Miskatonic](https://github.com/happy-Diode/Miskatonic) first, due to the absent of our internal network and authentication infrastructure, to protect QuEra packages, the current process of setting up the registry and internal packages is bit more complicated than the open source
workflow. Please send a message in the QuEra slack channel `#julia` if you get any questions.

## Add QuEra Private Registry

### Setting up GitHub authentication

In order for the Julia package manager to access this private registry, it will need to be able to
authenticate with Github.

Username-password authentication for git operations 
[has been deprecated](https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/):
- passwords are not accepted if two-factor auth is enabled for a user's account
- beginning on August 13, 2021 passwords will not be accepted for any git operations

Anyone who wishes to use this registry should generate a personal
access token and use it to add the registry to Julia:

1. configure git to [cache your credentials](https://docs.github.com/en/free-pro-team@latest/github/using-git/caching-your-github-credentials-in-git);
   refer to apporopriate OS-specific instructions

2. [generate a personal access token](https://github.com/settings/tokens); at a minimum, select the `repo`
   scope, and copy/paste the token into a scratch file

3. git clone a private repository from github with a HTTPS URL; when prompted to enter a password,
   paste the generated token

The last step will trigger caching of the token. The git library used by the Julia package manager
should be able to access this cache.

### Add the Miskatonic Registry

Please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode.

```jl
registry add https://github.com/Happy-Diode/Miskatonic.git
```

If you followed the Github authentication steps above to cache the personal access token,
this should add the registry without prompting for a username and password.

Please ask for help on Slack channel `#julia` if you encounter errors.
