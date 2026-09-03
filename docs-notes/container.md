# Container build and use

The root `Dockerfile` produces the same release `aie` binary as a local Cargo
build and runs it as an unprivileged user. Build it from the repository root:

```sh
docker build --pull -t gravlax:0.1.0 .
docker run --rm gravlax:0.1.0 --version
```

Mount a project at `/work`. Passing the host user and group lets plan outputs
remain host-writable on Linux:

```sh
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$PWD/analysis:/work" \
  gravlax:0.1.0 doctor --project /work

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$PWD/analysis:/work" \
  gravlax:0.1.0 plan run /work/plans/replay.yaml --project /work
```

The runtime image intentionally contains only `aie`. STAR and samtools are
optional doctor warnings and should be supplied in a separate alignment image
when needed; archive replay and queries do not require them.

For a read-only audit, add `--read-only` and mount the project read-only. The
workspace check will then report that outputs cannot be created, as intended.
