#!/bin/sh
set -ex
#RUSTFLAGS="-C target-cpu=native" cargo build --release --features print-trace --bin "expander-exec"
RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --bin expander-exec --release -- prove ./data/circuit_big_bn254.txt ./data/witness_big_bn254.txt ./data/proof.bin
RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --bin expander-exec --release -- verify ./data/circuit_big_bn254.txt ./data/witness_big_bn254.txt ./data/proof.bin
#RUSTFLAGS="-C target-cpu=native" cargo build --release --bin "expander-exec"
#RUST_LOG="warn" ./target/release/expander-exec prove ./data/circuit_small.txt ./data/witnessBinarySmall.txt ./data/proof.txt
#RUST_LOG="warn" ./target/release/expander-exec verify ./data/circuit_small.txt ./data/witnessBinarySmall.txt ./data/proof.txt

#RUST_LOG="warn" ./target/release/expander-exec prove ./data/circuit_big.txt ./data/witnessBinaryBig.txt ./data/proof.txt
#RUST_LOG="warn" ./target/release/expander-exec verify ./data/circuit_big.txt ./data/witnessBinaryBig.txt ./data/proof.txt
