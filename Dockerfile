# Ambuild's Poreblazer fork as a small runtime image: the executable at
# /opt/poreblazer/poreblazer.exe and the libraries it needs, on Debian bookworm.
# Published to ghcr.io/st7ma784/poreblazer by .github/workflows/ci.yml once the tests
# pass. Other images take the executable from it:
#
#   COPY --from=ghcr.io/st7ma784/poreblazer:sha-<commit> /opt/poreblazer/poreblazer.exe /opt/poreblazer/
#
# (they need libgfortran5 and libgomp1). Or run it on a directory holding input.dat,
# defaults.dat, UFF.atoms and the structure:
#
#   docker run --rm -i -v "$PWD:/work" ghcr.io/st7ma784/poreblazer:ambuild < input.dat

FROM debian:bookworm-slim AS build
RUN apt-get update \
 && apt-get install -y --no-install-recommends gfortran make \
 && rm -rf /var/lib/apt/lists/*
COPY src /src
# The Makefile's flags: -O2 -fopenmp -ffp-contract=off -fvect-cost-model=dynamic
RUN make -C /src && strip /src/poreblazer.exe

FROM debian:bookworm-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgfortran5 libgomp1 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/poreblazer.exe /opt/poreblazer/poreblazer.exe
COPY FORK.md LICENSE /opt/poreblazer/
LABEL org.opencontainers.image.source="https://github.com/st7ma784/poreblazer" \
      org.opencontainers.image.description="Ambuild's Poreblazer fork: OpenMP, cell lists, opt-in exact percolation labelling; output identical to Poreblazer 3.0.5 by default" \
      org.opencontainers.image.licenses="GPL-3.0"
ENV PATH=/opt/poreblazer:$PATH
WORKDIR /work
ENTRYPOINT ["poreblazer.exe"]
