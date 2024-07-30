#!/bin/sh
set -ex
RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --release --features print-trace --bin matmul 
#RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --release --bin matmul 
