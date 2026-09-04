# syntax=docker/dockerfile:1

FROM rust:1.89-bookworm AS builder
WORKDIR /src
COPY . .
RUN cargo build --locked --release -p gravlax

FROM debian:bookworm-slim AS runtime
LABEL org.opencontainers.image.title="Gravlax" \
      org.opencontainers.image.description="Annotation-independent molecular evidence archive and query CLI" \
      org.opencontainers.image.source="https://github.com/COMBINE-lab/gravlax" \
      org.opencontainers.image.licenses="BSD-3-Clause"

COPY --from=builder /src/target/release/aie /usr/local/bin/aie
RUN mkdir /work && chown 65532:65532 /work

USER 65532:65532
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/aie"]
CMD ["--help"]
