#!/bin/sh
set -ex

if [ $1 ]
then
        echo "Debug mode"
        # debug
        RUST_LOG="debug" cargo run --bin matmul
else
        echo "Release mode"
        # release
        #RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --release --features print-trace --bin matmul 
        RUST_LOG="info" RUSTFLAGS="-C target-cpu=native" cargo run --release --bin matmul 
fi