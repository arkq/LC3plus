FROM --platform=linux/amd64 container-registry.gitlab.cc-asp.fraunhofer.de/lc3-lc3plus-extern/lc3_etsi_release/lc3plus_conformance:1.0.18

# Fix Debian Buster EOL repository issue - switch to archive
RUN sed -i 's|http://deb.debian.org/debian|http://archive.debian.org/debian|g' /etc/apt/sources.list && \
    sed -i 's|http://security.debian.org/debian-security|http://archive.debian.org/debian-security|g' /etc/apt/sources.list && \
    sed -i '/buster-updates/d' /etc/apt/sources.list

# Install 32-bit MinGW cross-compiler for Windows
RUN apt-get update && apt-get install -y \
    gcc-mingw-w64-i686 \
    && rm -rf /var/lib/apt/lists/*

# Add existing Rust to PATH
ENV PATH="/root/.cargo/bin:${PATH}"

# Add 32-bit Windows target
RUN rustup target add i686-pc-windows-gnu

# Install Python packages
RUN pip install scipy soundfile

WORKDIR /app