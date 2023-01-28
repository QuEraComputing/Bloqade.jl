release-patch:
    #!/usr/bin/env bash
    set -e
    git pull origin master
    git push origin master

    for i in lib/*; do
        if [ -d "$i" ]; then
            cd $i
            ion bump patch --no-commit
            git add Project.toml
            cd ../..
        fi
    done
    ion bump patch --no-commit
    git commit -m "Bump patch version"
    git push origin master

    for i in lib/*; do
        if [ -d "$i" ]; then
            ion summon lib/$i --skip-note
        fi
    done
    ion summon

release-minor:
    error("Please manually bump the minor version using ion and change corresponding compat")
