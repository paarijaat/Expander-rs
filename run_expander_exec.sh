#!/bin/sh
set -ex
#RUSTFLAGS="-C target-cpu=native" cargo build --release --features print-trace --bin "expander-exec"
RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --release --bin expander-exec -- prove ./data/circuit_big_bn254.txt ./data/witness_big_bn254.txt ./data/proof.bin
RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --release --bin expander-exec -- verify ./data/circuit_big_bn254.txt ./data/witness_big_bn254.txt ./data/proof.bin
