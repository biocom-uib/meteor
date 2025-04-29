podman run -v "$PWD:/volume" --rm -t -w /volume docker.io/clux/muslrust:1.88.0-nightly cargo build --release
