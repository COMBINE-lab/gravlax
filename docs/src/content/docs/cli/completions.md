---
title: Shell completions
description: Generate current command and option completions for Bash, Zsh, or Fish.
---

`aie completions` generates a script from the command interface in the
installed binary. This keeps completions aligned with the exact Gravlax
version in use.

## Bash

For the current shell:

```sh
source <(aie completions bash)
```

For future sessions, write the generated script into a directory loaded by
your Bash completion setup, commonly:

```sh
mkdir -p ~/.local/share/bash-completion/completions
aie completions bash > ~/.local/share/bash-completion/completions/aie
```

## Zsh

Choose a directory on `fpath`, then generate `_aie` there:

```sh
mkdir -p ~/.zfunc
aie completions zsh > ~/.zfunc/_aie
```

Add `~/.zfunc` to `fpath` before `compinit` in `.zshrc` if it is not already
configured.

## Fish

```sh
mkdir -p ~/.config/fish/completions
aie completions fish > ~/.config/fish/completions/aie.fish
```

Regenerate the file after upgrading Gravlax so new commands and options are
included.
