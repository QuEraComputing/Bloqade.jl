release-patch:
    ion bump patch lib/BloqadeExpr            --no-commit
    ion bump patch lib/BloqadeKrylov          --no-commit
    ion bump patch lib/BloqadeLattices        --no-commit
    ion bump patch lib/BloqadeMIS             --no-commit
    ion bump patch lib/BloqadeODE             --no-commit
    ion bump patch lib/BloqadeSchema          --no-commit
    ion bump patch lib/BloqadeWaveforms       --no-commit
    ion bump patch lib/YaoSubspaceArrayReg    --no-commit
    ion bump patch

    git add lib/BloqadeExpr/Project.toml
    git add lib/BloqadeKrylov/Project.toml
    git add lib/BloqadeLattices/Project.toml
    git add lib/BloqadeMIS/Project.toml
    git add lib/BloqadeODE/Project.toml
    git add lib/BloqadeSchema/Project.toml
    git add lib/BloqadeWaveforms/Project.toml
    git add lib/YaoSubspaceArrayReg/Project.toml
    git commit -m "Bump patch version"

    ion summon

release-minor:
    error("Please manually bump the minor version using ion and change corresponding compat")
